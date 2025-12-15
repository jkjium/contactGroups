import os.path
import time
import random
import itertools

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import typing as T

import torch
import torch.utils.data as data
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torchmetrics as tm
import torch.optim as optim
import torch.nn.functional as F
from torch.nn.utils.rnn import pad_sequence

import transformers, datasets
from transformers import EsmModel, AutoTokenizer

import peft
from peft import LoraConfig, TaskType, get_peft_model

transformers.logging.set_verbosity_error()

from pathlib import Path
from tqdm import tqdm
import commp as cp
import logging

ESMs = [ "facebook/esm2_t6_8M_UR50D" ,
         "facebook/esm2_t12_35M_UR50D" ,
         "facebook/esm2_t30_150M_UR50D" ,
         "facebook/esm2_t33_650M_UR50D" ,
         "facebook/esm2_t36_3B_UR50D" ]

logfile = f"{Path(__name__).resolve().parent.name}_.log"
logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s | %(levelname)s | %(message)s',
    handlers=[
        logging.FileHandler(logfile),
        logging.StreamHandler()               # Also print to console
    ]
)
logger = logging.getLogger(__name__)
logger.info(f"Logging to: {logfile}")

##############################################################
# ESM + Lora + cls_head
class cls_head(nn.Module):
    """Head for sentence-level classification tasks."""

    def __init__(self, input_size, hidden_size, hidden_dropout_prob, n_classes):
        super().__init__()
        self.dense = nn.Linear(input_size, hidden_size)
        self.dropout = nn.Dropout(hidden_dropout_prob)
        self.out_proj = nn.Linear(hidden_size, n_classes)

    def forward(self, features, **kwargs):
        x = self.dropout(features)
        x = self.dense(x)
        x = torch.tanh(x)
        x = self.dropout(x)
        x = self.out_proj(x)
        return x

class cls_model(nn.Module):
    def __init__(self, language_model, classifier):
        super().__init__()
        self.language_model = language_model
        self.classifier = classifier
        
    def forward(self, input_ids, attention_mask):
        _, lm_embedding = self.language_model(input_ids, attention_mask)
        clsf = self.classifier(lm_embedding)
        return clsf


class PPIModel(nn.Module):
    def __init__(self, language_model, classifier):
        super().__init__()
        self.language_model = language_model
        self.classifier = classifier
        
    def forward(self, input_ids, attention_mask):
        p1_tokens, p2_tokens = input_ids[0], input_ids[1]
        p1_attn, p2_attn = attention_mask[0], attention_mask[1]
        a = self.language_model(p1_tokens, p1_attn)[1]
        b = self.language_model(p2_tokens, p2_attn)[1]
        order_invariant_emb = torch.cat((torch.abs(a - b), a * b), dim=1)
        return self.classifier(order_invariant_emb)
        #mean_embedding = torch.div(torch.add(p1_embedding, p2_embedding), 2)
        #return self.classifier(mean_embedding)
        

class PPIClassificationHead(nn.Module):
    def __init__(self, input_size, hidden_size, hidden_dropout_prob):
        super().__init__()
        self.dense = nn.Linear(input_size, hidden_size)
        self.dropout = nn.Dropout(hidden_dropout_prob)
        self.out_proj = nn.Linear(hidden_size, 1)

    def forward(self, features, **kwargs):
        x = self.dropout(features)
        x = self.dense(x)
        x = torch.tanh(x)
        x = self.dropout(x)
        x = self.out_proj(x)
        return x


# functions
def add_weight_decay(
    model: nn.Module,
    l2_coeff: float,
    lora_l2_coeff: T.Optional[float] = None
    ) -> T.List[T.Dict[str, T.Union[nn.Parameter, float]]]:

    lora_l2_coeff = lora_l2_coeff or l2_coeff

    lora_decay, decay, no_decay = [], [], []
    for name, param in model.named_parameters():
        if not param.requires_grad:
            continue
        elif "lora" in name:
            lora_decay.append(param)
        elif "norm" in name or name.endswith(".bias"):
            no_decay.append(param)
        else:
            decay.append(param)
    return [{'params': lora_decay, 'weight_decay': lora_l2_coeff}, {'params': no_decay, 'weight_decay': 0.0}, {'params': decay, 'weight_decay': l2_coeff}]


## do_peft parameters: 
## layers_to_transform:
##   - 32
##   - 31
##   - 30
##   - 29
##   - 28
##   - 27
##   - 26
##   - 25
## inference_mode: false
## lora_r: 8
## lora_alpha: 32
## lora_dropout: 0.1
## lora_bias: none
## lora_l2_coeff: 0.01
## do_query: False
## do_key: True
## do_value: True

# peft ft procedure - get_model
## set peft parameters
## load esm model
## generate peft injected model
## add classifier head
## setup trainable grads for whole model
## return model
def esm_lora_cls_model():
    ## https://huggingface.co/docs/peft/main/en/package_reference/lora#peft.LoraConfig.target_modules
    peft_config = LoraConfig(
        task_type=TaskType.FEATURE_EXTRACTION, ## extracting embeddings or representations from the model (not classification or )
        inference_mode=False, ## {True: LoRA weights are frozen, and only used during inference (no training), False: LoRA layers are trainable (used for fine-tuning)}
        r=8, ## Lora rank, full size d x d, then Lora using two matrices: d x r and r x d; smaller r means fewer trainable parameters
        lora_alpha=32, ## scaling factor = alpha / r; controls how much Lora weights affect the base model; higher value means stronger influence
        lora_dropout=0.1, ## (only applies during training) 0.0 for inference, 0.05 - 0.1 for training
        bias="none", ## usually none {all, lora_only}
        target_modules=['key', 'value'], ## The names of the modules to apply the adapter to. print(model) will show names. 
        ## ESM has 
        ##    (query): Linear(in_features=1280, out_features=1280, bias=True)
        ##    (key): Linear(in_features=1280, out_features=1280, bias=True)
        ##    (value): Linear(in_features=1280, out_features=1280, bias=True)
        ##    (dropout): Dropout(p=0.0, inplace=False)
        ##    (rotary_embeddings): RotaryEmbedding()
        layers_to_transform=[32, 31, 30, 29, 28, 27, 26, 25], ## (Optional) List of layer indices to apply LoRA to. Reduces computation and focuses learning on impactful layers
        modules_to_save=[], ## module names in base model (other than LoRA adapter modules) to keep trainable and save
    )

    plm = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D", return_dict = False)
    base_model = get_peft_model(plm, peft_config)

    ## In Hugging Face's EsmModel, the pooler is used to summarize the entire sequence into a fixed-size vector, usually based on the [CLS] token
    ## 1. take [CLS] token (size of 1280) from the last hidden layer; 2. pass to 1280 x 1280 linear layer and activation function; output size 1280
    input_size = plm.pooler.dense.out_features * 2 # 2560 1280 * 2
    print(input_size)
    ## plm.pooler = EsmPooler(config)
    ## EsmPooler:
    ##      self.dense = nn.Linear(config.hidden_size, config.hidden_size)
    ##      self.activation = nn.Tanh() 
    #classifier_head = multiClassificationHead(input_size = input_size, hidden_size = 128, hidden_dropout_prob = 0.0, n_classes = 18) # from config
    classifier_head = PPIClassificationHead(input_size=input_size, hidden_size=128, hidden_dropout_prob = 0.0) # from config
    model = PPIModel(language_model = base_model, classifier = classifier_head)
    model = model.to('cuda', non_blocking = True)    

    # set parameters to be trained
    for (name, param) in model.named_parameters():
        param.requires_grad = False
        for k in ['classifier', 'lora']:
            if k in name:
                param.requires_grad = True
                break

    #if inference_mode:
    #    model = model.eval()
    return model


# data loading procedure
# tokenization within dataset 
# easier for fine-tuning
class ESMTokenizerDataset(Dataset):
    def __init__(self, data_df: pd.DataFrame, max_crop: int = 1024, esm_variant: str = "facebook/esm2_t33_650M_UR50D"):
        self.data_df = data_df
        self.max_seq_len = max_crop - 2
        self.tokenizer = AutoTokenizer.from_pretrained(esm_variant)
    
    def __len__(self):
        return len(self.data_df) # number of sequence pairs

    def __getitem__(self, index: int) -> T.Tuple[torch.Tensor, torch.Tensor]:
        row = self.data_df.iloc[index] # s0, s1, label
        s0, s1 = row[0][:self.max_seq_len], row[1][:self.max_seq_len] # crop sequences > max_length
        tok0 = self.tokenizer(s0, return_tensors="pt")["input_ids"].squeeze()
        tok1 = self.tokenizer(s1, return_tensors="pt")["input_ids"].squeeze()
        label = torch.tensor([int(row[2])]) # interact or not
        return (tok0, tok1, label)

# collate function for padding token
# for batch_size = 4:
#   batch len(b) = 3 
#   b[0][0] = torch.Size([4, 1024]) # seq tok0
#   b[0][1] = torch.Size([4, 1024]) # seq tok1
#   b[1] = torch.Size([4, 1]) # label
#   b[2][0] = torch.Size([4, 1024]) # seq tok mask0
#   b[2][1] = torch.Size([4, 1024]) # seq tok mask1
def _fn_loader_collate(batch):
    tok0, tok1, labels = zip(*batch)
    # PAD tokens with PAD_VALUE to have the same length
    # batch_first=True: Output tensor shape:[batch_size, max_seq_len, feature_dim] instead of [max_seq_len, batch_size, feature_dim]
    tok0 = pad_sequence(tok0, batch_first = True, padding_value = 0.)
    tok1 = pad_sequence(tok1, batch_first = True, padding_value = 0.)
    labels = torch.stack(labels)
    attn_mask0 = (tok0!=0).int()
    attn_mask1 = (tok1!=0).int()
    return ((tok0, tok1), labels, (attn_mask0, attn_mask1)) 

# dataset instance + data loader setup
# max_crop = 1024
# batch_size = 4
def esm_lora_cls_load_data(params):
    # torch Dataset encapsulation 
    train_df = pd.read_csv('01.train.seq2.label.tsv', header=None, sep='\t', dtype='str')
    train_dataset = ESMTokenizerDataset(train_df, 1024, "facebook/esm2_t33_650M_UR50D")
    valid_df = pd.read_csv('01.valid.seq2.label.tsv', header=None, sep='\t', dtype='str')
    valid_dataset = ESMTokenizerDataset(valid_df, 1024, "facebook/esm2_t33_650M_UR50D")

    # torch Data loader
    # pin_memory: faster data loading when training; num_workers: num of CPUs parallel load data
    train_loader = DataLoader(train_dataset, batch_size=params.batch_size, shuffle=True, num_workers=0, pin_memory=True, collate_fn=_fn_loader_collate)
    valid_loader = DataLoader(valid_dataset, batch_size=params.batch_size, shuffle=False, num_workers=0, pin_memory=True, collate_fn=_fn_loader_collate)    
    return train_loader, valid_loader


# metric collection factory function
# with each batch: m.reset() -> m.update(preds, targets) -> m.compute()
def get_metric_collection(prefix: str = "valid/"):
    collection = tm.MetricCollection({
        "accuracy": tm.Accuracy(task = "binary"),
        "aupr": tm.AveragePrecision(task = "binary"),
        "precision": tm.Precision(task = "binary"),
        "recall": tm.Recall(task = "binary"),
        "f1": tm.F1Score(task = "binary"),
        "mcc": tm.MatthewsCorrCoef(task = "binary"),
        "specificity": tm.Specificity(task = "binary"),
    }, prefix = prefix)
    return collection


def save_model(checkpoint_key, model, optimizer, scheduler, epoch, best_valid_loss, best_valid_aupr, params):
    ## distributed version 
    #state_dict = model.module.state_dict()
    state_dict = model.state_dict()
    
    for k,v in state_dict.items():
        if isinstance(v, torch.Tensor):
            state_dict[k] = v.cpu()
    save_dict = {
        "epoch": epoch,
        "model_state_dict": state_dict,
        "optimizer_state_dict:": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "best_valid_loss": best_valid_loss, 
        "best_valid_aupr": best_valid_aupr, 
    }
    model_path = f"./models/{params.run_name}"
    Path(model_path).mkdir(parents=True, exist_ok=True)
    checkpoint_name = f"{model_path}/{params.run_name}_{params.esm_pretrained}_finetune_{checkpoint_key}.pt"
    torch.save(save_dict, checkpoint_name)

# procedure:
# load data
# create model
# setup optimizer
# setup scheduler
# train
# parameters
# lora_l2_coeff: 0.01, l2_coeff: 0.01, lr: 0.001
def esm_lora_cls_train(params):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f'batch size: {params.batch_size}')

    # load data
    train_loader, valid_loader = esm_lora_cls_load_data(params)
    # create model
    model = esm_lora_cls_model()
    print(model)
    model.to(device)

    # init torch metrics
    tr_metrics = get_metric_collection(prefix = "train/").to(device)
    tr_metrics.reset()
    vl_metrics = get_metric_collection(prefix = "valid/").to(device)
    vl_metrics.reset()
    active_fn = nn.Sigmoid() # for metrics of logits

    # set optimizer
    # (optional) add weight decay
    # using different decay for different parameters
    # opt_params: return [{'params': lora_decay, 'weight_decay': lora_l2_coeff}, {'params': no_decay, 'weight_decay': 0.0}, {'params': decay, 'weight_decay': l2_coeff}]
    opt_params = add_weight_decay(model, params.l2_coeff, params.lora_l2_coeff)
    # non decay optimizer
    # optimizer = torch.optim.AdamW(model.parameters(), lr=params.lr, eps=params.epsilon)
    optimizer = optim.AdamW(opt_params, lr=params.lr) # torch.optim as optim
    # optional
    # T_0=10 first cycle length 10 epochs; T_mult=1, how cycle length changes over time: T_n = T_mult * T_(n-1), n>1; eta_min=4e-6, minimum LR value
    scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, 10, 1, 4e-6)
    loss_fn = nn.BCEWithLogitsLoss() # for binary prediction
    logger.info('model setup done.')

    # load checkpoint 
    if hasattr(params, "checkpoint"):
        logger.info(f"Loading checkpoint: {params.checkpoint}")
        
        checkpoint = torch.load(params.checkpoint, map_location = device) #torch.device(params.rank)
        model.module.load_state_dict(checkpoint["model_state_dict"])
        optimizer.load_state_dict(checkpoint["optimizer_state_dict"])
        scheduler.load_state_dict(checkpoint["scheduler_state_dict"])
        start_epoch = checkpoint["epoch"] + 1

        best_valid_loss = checkpoint["best_valid_loss"]
        best_valid_aupr = checkpoint["best_valid_aupr"]
    else:
        start_epoch = 0
        best_valid_loss = torch.inf
        best_valid_aupr = 0

    # epoch loop
    tr_loss_log = 0
    vl_loss_log = 0
    for epoch in range(start_epoch, params.n_epoch):
        # training loop
        model.train()
        optimizer.zero_grad() # redundant, defensive coding
        for i, batch in enumerate(train_loader):
            #logger.info(f'epoch: {epoch} trianing step: {i} / {len(train_loader)}')
            with torch.amp.autocast('cuda', enabled=params.USE_AMP): # autocast precision, improve speed
                (seq_tok0, seq_tok1), labels, (attn0, attn1) = batch
                feats = (seq_tok0.to(device), seq_tok1.to(device))
                labels = labels.to(device)
                attn = (attn0.to(device), attn1.to(device))                

                logits = model(feats, attention_mask=attn)
                loss = loss_fn(logits, labels.float()) / params.accum_step # normalize loss to avoid gradient explosion
                tr_loss_log += loss.detach()
                loss.backward()

                tr_metrics(active_fn(logits).float().detach(), labels.int()) # torch.metrics.update()

            if (i+1) % params.accum_step == 0 or (i+1) == len(train_loader): # gradiant accumulation; to achieve big batch size on a small GPU
                torch.nn.utils.clip_grad_norm_(model.parameters(), 0.2)
                optimizer.step()
                optimizer.zero_grad()
                scheduler.step() # update lr in optimizer
                logger.info(f"Current lr: {optimizer.param_groups[0]['lr']}")

            if ((i + 1) % params.report_step == 0):
                report_loss = tr_loss_log.item() * params.accum_step / (i+1)
                logger.info(f"[Train] [Epoch {epoch}] Step {i+1}/{len(train_loader)} | Loss: {report_loss:.6f}")


        # output metrics
        train_ret = tr_metrics.compute()
        train_ret['tr_loss'] = tr_loss_log
        for k, v in train_ret.items():
            logger.info(f"[Metric] [Step {epoch}] {k} = {v:.6f}")

        # validating loop
        model.eval() # freeze dropout and BatchNorm layers
        for i, batch in enumerate(valid_loader):
            #logger.info(f'epoch: {epoch} validating step: {i} / {len(valid_loader)}')
            with torch.no_grad(): # no training
                with torch.amp.autocast('cuda', enabled=params.USE_AMP):
                    (seq_tok0, seq_tok1), labels, (attn0, attn1) = batch
                    feats = (seq_tok0.to(device), seq_tok1.to(device))
                    labels = labels.to(device)
                    attn = (attn0.to(device), attn1.to(device))

                    logits = model(feats, attention_mask=attn)
                    loss = loss_fn(logits, labels.float())
                    vl_loss_log += loss.detach()

                    vl_metrics(active_fn(logits).float().detach(), labels.int()) # same as metrics.update(predicts, targets)

                if ((i + 1) % params.report_step == 0):
                    logger.info(f"[Validation] [Epoch {epoch}] Step {i+1}/{len(valid_loader)}")                

        # output metrics
        valid_ret = vl_metrics.compute()
        valid_ret['valid/loss'] = vl_loss_log / len(valid_loader)
        for k, v in valid_ret.items():
            logger.info(f"[Metric] [Step {epoch}] {k} = {v:.6f}")

        # procedure save model
        if valid_ret["valid/loss"] < best_valid_loss:
            logger.info(f"Save model: New best validation loss: {valid_ret['valid/loss']:.6f}")
            best_valid_loss = valid_ret["valid/loss"]
            save_model("best_loss", model, optimizer, scheduler, epoch, best_valid_loss, best_valid_aupr, params)

        if valid_ret["valid/aupr"] > best_valid_aupr:
            logger.info(f"Save model: New best validation AUPR: {valid_ret['valid/aupr']:.6f}")
            best_valid_aupr = valid_ret["valid/aupr"]
            save_model("best_aupr", model, optimizer, scheduler, epoch, best_valid_loss, best_valid_aupr, params)
            
        if (epoch + 1) % params.save_every == 0:
            logger.info(f"Saving checkpoint (Epoch {epoch})")
            save_model(f"checkpoint{epoch}", model, optimizer, scheduler, epoch, best_valid_loss, best_valid_aupr, params)

    save_model("last", model, optimizer, scheduler, epoch, best_valid_loss, best_valid_aupr, params)



def esm_lora_cls_main():
    from omegaconf import OmegaConf
    print(torch.cuda.device_count())
    params = OmegaConf.load('config/ppi.yaml')
    os.makedirs(f"./models/{params.run_name}", exist_ok=True)
'''
import os.path
import time
import random
import itertools
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import typing as T

import torch
import torch.utils.data as data
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torchmetrics as tm
import torch.nn.functional as F
from torch.nn.utils.rnn import pad_sequence
import transformers, datasets
from transformers import EsmModel, AutoTokenizer
import peft
from peft import LoraConfig, TaskType, get_peft_model
transformers.logging.set_verbosity_error()

from tqdm import tqdm
from pathlib import Path
import commp as cp
import importlib as ib
import utils_plm as up
from omegaconf import OmegaConf

#logger = logging.getLogger(__name__)
#logger.info(f"Logging to: {logfile}")
params = OmegaConf.load('config/ppi.yaml')

up.esm_lora_cls_train(params)


'''
##############################################################
# multi-fc layer
# Input → {Linear → ReLU → Dropout} x n → Output
class predict_layer(nn.Module):
    def __init__(self, input_dim, dense_units, dropout_rate, num_label):
        super(predict_layer, self).__init__()
        self.normalizer = nn.BatchNorm1d(input_dim)
        self.fc1 = nn.Linear(input_dim, dense_units)
        self.dropout = nn.Dropout(dropout_rate)
        self.output = nn.Linear(dense_units, num_label)
        self.num_label = num_label

    def forward(self, x):
        x = self.normalizer(x)
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.output(x)  # raw output (logits)
        return x

class predict_layer2(nn.Module):
    def __init__(self, input_dim, dense_units, dropout_rate, num_label):
        super(predict_layer2, self).__init__()
        self.normalizer = nn.BatchNorm1d(input_dim)
        self.fc1 = nn.Linear(input_dim, dense_units)
        self.dropout1 = nn.Dropout(dropout_rate)

        # self.normalizer2 = nn.BatchNorm1d(dense_units) 
        # previous self.normalizer is tuned for the input dimension only
        # adding another normalizier does not increase performance
        self.fc2 = nn.Linear(dense_units, dense_units)
        self.dropout2 = nn.Dropout(dropout_rate)

        self.output = nn.Linear(dense_units, num_label)
        self.num_label = num_label

    def forward(self, x):
        x = self.normalizer(x)
        x = F.relu(self.fc1(x))
        x = self.dropout1(x)

        # x = self.normalizer2(x)
        x = F.relu(self.fc2(x))
        x = self.dropout2(x)

        x = self.output(x)  # raw output (logits)
        return x

# generalized multi-layer predictor
class mpredict_layer(nn.Module):
    def __init__(self, input_dim, dense_units, dropout_rate, num_label, num_hidden_layers=1):
        super(mpredict_layer, self).__init__()

        self.normalizer = nn.BatchNorm1d(input_dim)

        layers = []
        in_dim = input_dim
        for _ in range(num_hidden_layers):
            layers.append(nn.Linear(in_dim, dense_units))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout_rate))
            in_dim = dense_units  # output of this layer is input to next

        self.hidden_layers = nn.Sequential(*layers)
        self.output = nn.Linear(dense_units, num_label)
        self.num_label = num_label

    def forward(self, x):
        x = self.normalizer(x)
        x = self.hidden_layers(x)
        x = self.output(x)
        return x

class TabularDataset(data.Dataset):
    def __init__(self, df, target_col):
        self.X = torch.tensor(df.iloc[:, :-2].values, dtype=torch.float32)
        self.y = torch.tensor(df.iloc[:, target_col].values, dtype=torch.float32)
        if len(self.y.shape) == 1:
            self.y = self.y.unsqueeze(1)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

###################################################
# w/o ft procedure
# input data: protein embeddings and labels
###################################################
# 
# spearman, preds = up.test_predictor(m)
# spearman: 0.8030012634744585 
def test_predictor(model):
    df = pd.read_pickle("./notebooks/embedding/example_data/GB1_embedded/test_ESM2_8M.pkl")
    target_col = -1
    model.eval()
    x_test = torch.tensor(df.iloc[:, :-2].values, dtype=torch.float32).to(next(model.parameters()).device)
    y_test = df.iloc[:, target_col].values
    # next(model.parameters()).device::Gets the first parameter tensor to determine the device of parameters  (one is enough)
	
    with torch.no_grad():
        preds = model(x_test).cpu().numpy().squeeze()
		
    spearman = spearmanr(y_test, preds).correlation
    return spearman, preds

# model, spearman_scores = up.train_predictor(seed = 42)    
def train_predictor(epochs=240, lr=1e-4, epsilon=1e-7, batch=8, dropout=0.2, dense=32, seed=99, num_labels=1, num_hidden_layers = 1):
    # Set seed
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

    # Load and split your data (placeholder - replace `read()` with your own function)
    df_train=pd.read_pickle("./notebooks/embedding/example_data/GB1_embedded/train_ESM2_8M.pkl")
    df_valid=pd.read_pickle("./notebooks/embedding/example_data/GB1_embedded/valid_ESM2_8M.pkl")
    df_test=pd.read_pickle("./notebooks/embedding/example_data/GB1_embedded/test_ESM2_8M.pkl")

    # Create datasets and loaders
    train_dataset = TabularDataset(df_train, target_col=-1)
    valid_dataset = TabularDataset(df_valid, target_col=-1)

    train_loader = data.DataLoader(train_dataset, batch_size=batch, shuffle=True)
    valid_loader = data.DataLoader(valid_dataset, batch_size=batch, shuffle=False)

    # Model initialization
    input_dim = df_train.shape[1] - 2  # Assuming last 2 columns are labels/metadata
    print('num_hidden_layer: %d' % num_hidden_layers)
    #model = predict_layer(input_dim=input_dim, dense_units=dense, dropout_rate=dropout, num_label=num_labels)
    #model = predict_layer2(input_dim=input_dim, dense_units=dense, dropout_rate=dropout, num_label=num_labels)
    model = mpredict_layer(input_dim=input_dim, dense_units=dense, dropout_rate=dropout, num_label=num_labels, num_hidden_layers = num_hidden_layers)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)

    # Optimizer and loss
    print('Using AdamW.')
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, eps=epsilon)
    loss_fn = nn.MSELoss() if num_labels == 1 else nn.CrossEntropyLoss()

    # Track spearman scores
    spearman_scores = []

    # Training loop
    for epoch in range(epochs):
        model.train()
        for xb, yb in train_loader:
            xb, yb = xb.to(device), yb.to(device)

            optimizer.zero_grad()
            preds = model(xb)

            if num_labels == 1:
                loss = loss_fn(preds, yb)
            else: # categorical 
                loss = loss_fn(preds, yb.long().squeeze())

            loss.backward()
            optimizer.step()

        # Validation Spearman calculation
        model.eval()
        all_preds = []
        all_targets = []

        with torch.no_grad():
            for xb, yb in valid_loader:
                xb = xb.to(device)
                preds = model(xb).cpu().numpy()
                all_preds.append(preds)
                all_targets.append(yb.numpy())

        all_preds = np.vstack(all_preds).squeeze()
        all_targets = np.vstack(all_targets).squeeze()
        
        # Compute Spearman correlation
        score = spearmanr(all_targets, all_preds).correlation
        spearman_scores.append(score)

        print(f"Epoch {epoch+1}/{epochs} - Spearman: {score:.4f}")

    return model, spearman_scores


# df_seq: {columns: sequence, ....}
def _emb_esm(df, emb_type='prot', checkpoint='facebook/esm2_t33_650M_UR50D'):
    # setup model
    tokenizer = AutoTokenizer.from_pretrained(checkpoint)
    model = EsmModel.from_pretrained(checkpoint, torch_dtype=torch.float16)
    model = model.to("cuda")
    model = model.half()
    
    # embedding
    emb = []
    for i in tqdm(range(0,len(df))):
        inputs = tokenizer(df["sequence"].loc[i], return_tensors="pt", max_length = 10000, truncation=True, padding=False).to("cuda")
        with torch.no_grad():
            if emb_type == 'prot':
                # .cpu() required for converting to NumPy
                emb.append(np.array(torch.mean( model(**inputs).last_hidden_state.cpu(), dim = 1)))
            else:
                out = np.array( model(**inputs).last_hidden_state.cpu()) # out = [batch_size, seq_len, hidden_dim]
                out = np.squeeze(out) # remove singleton dimension. here is to remove batch_size dimention if batch_size = 1
                out = out[1:-1, :] # each embedding: [CLS, res_emb, SEP], out[1:-1,:] remove the first and last special token embeddings
                emb.append(out)
    return emb
    #return df_emb

def test_emb(infile = './training data/SecStr/test.pkl', emb_type='prot'):
    df_seq = pd.read_pickle(infile)
    df_seq["sequence"]=df_seq["sequence"].str.replace('|'.join(["O","B","U","Z","J"]),"X",regex=True)
    df = df_seq[['sequence']].head(2)
    print(len(df['sequence'][0]), df['sequence'][0])
    print(len(df['sequence'][1]), df['sequence'][1])

    ret = _emb_esm(df, emb_type=emb_type)
    print('\n'.join(['%s: %s' % (str(i), str(ret[i].shape)) for i in range(len(ret))]))
    return ret
    #df_emb = pd.DataFrame(np.concatenate(ret))
    # concatenating multiple sources can lead to non-sequential indices
    # reset_index set a new sequential index; drop =True: remove the old index
    #df_emb.reset_index(drop = True, inplace = True)


if __name__=='__main__':
    cp.dispatch(__name__)