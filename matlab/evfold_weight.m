function W = evfold_weight(msa_fasta_filename, seqid_of_interest, theta)
[encoded_focus_alignment, focus_index_of_interest, focus_to_uniprot_offset_map] = read_alignment_fasta(msa_fasta_filename, seqid_of_interest);
[alignment_height,alignment_width] = size(encoded_focus_alignment);
W = ones(1, alignment_height);
if(theta > 0.0)
	W = (1./(1+sum(squareform(pdist(encoded_focus_alignment, 'hamm')<theta))));
end
end


function [encoded_focus_alignment, focus_index_of_interest, focus_to_uniprot_offset_map] = read_alignment_fasta(msa_fasta_filename, seqid_of_interest)
METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = getenv('DI_METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES'); % 1 = change them to gaps .. 2 = mask entire sequence
%if (size(METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES,2) == 0)
%	METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = 2;
%else
%	METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = str2num(METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES)
%end
METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES = 1;
full_alignment = fastaread(msa_fasta_filename);
alignment_width = size(full_alignment(1).Sequence, 2);
alignment_height = size(full_alignment, 1);
letter2number_map = create_letter2number_map();
[full_index_of_interest, range_of_interest_start, range_of_interest_end] = find_seq_of_interest(full_alignment, seqid_of_interest);
encoded_focus_alignment = [];
skipped_sequence_counter = 0;
[focuscolumnlist, focus_to_uniprot_offset_map] = scan_sequence_of_interest_for_focus_columns(full_alignment(full_index_of_interest).Sequence, range_of_interest_start, letter2number_map);
for full_alignment_index=1:alignment_height
	focus_alignment_row = full_alignment(full_alignment_index).Sequence(focuscolumnlist);
	encoded_focus_alignment_row = letter2number_map(focus_alignment_row);
	if (size(find(encoded_focus_alignment_row == 0),2) > 0)
		error(['Error: sequence in alignment has illegal characters: ' full_alignment(full_alignment_index).Sequence]);
	end
	%if (size(find(encoded_focus_alignment_row <= -2),2) > 0)
	%	error(['Error: sequence in alignment has dot or lowercase in conserved position: ' full_alignment(full_alignment_index).Sequence]);
	%end
	if (size(find(encoded_focus_alignment_row == -1),2) > 0)
		if (METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES == 1)
			encoded_focus_alignment_row(find(encoded_focus_alignment_row == -1)) = 1;
		else
			if (METHOD_TO_RESOLVE_AMBIGUOUS_RESIDUES == 2)
				continue %skip sequences with ambiguous residues
			else
				error('Internal Error');
			end
		end
	end
	encoded_focus_alignment(size(encoded_focus_alignment,1) + 1,:) = encoded_focus_alignment_row;
	if (full_alignment_index == full_index_of_interest)
		focus_index_of_interest = size(encoded_focus_alignment,1);
	end
end
end

function [index_of_interest, range_start, range_end] = find_seq_of_interest(full_alignment, seqid_of_interest)
index_of_interest = -1;
for scan_index = 1:size(full_alignment,1)
	[seqid, range_start, range_end] = split_uniprot_id(full_alignment(scan_index).Header);
	if (strcmp(seqid,seqid_of_interest) == 1)
		index_of_interest = scan_index;
		break
	end
end
if (index_of_interest == -1)
	error(['Error: could not find sequence of interest (' seqid_of_interest ') in multiple sequence alignment']);
end
end


function [seqid, range_start, range_end] = split_uniprot_id(pfam_uniprot_range_line)
slashposition = findstr('/', pfam_uniprot_range_line);
if (size(slashposition,2) ~= 1 || slashposition == 1 || slashposition == size(pfam_uniprot_range_line,2))
	error(['Error: could not parse (slash error) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
seqid = pfam_uniprot_range_line(1:slashposition - 1);
rangestring = pfam_uniprot_range_line(slashposition + 1:size(pfam_uniprot_range_line,2));
hyphenposition = findstr('-', rangestring);
if (size(hyphenposition,2) ~= 1 || hyphenposition == 1 || hyphenposition == size(rangestring,2))
	error(['Error: could not parse (hyphen error) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
range_start = str2num(rangestring(1:hyphenposition - 1));
range_end = str2num(rangestring(hyphenposition + 1 : size(rangestring,2)));
if (isempty(range_start) || isempty(range_end))
	error(['Error: could not parse (range start/end) uniprot range line in pfam alignment : ' pfam_uniprot_range_line]);
end
end


function [focuscolumnlist, uniprotoffsetlist] = scan_sequence_of_interest_for_focus_columns(sequence_of_interest, range_of_interest_start, letter2number_map)
focuscolumnlist = [];
uniprotoffsetlist = [];
next_uniprotoffset = range_of_interest_start;
for pos=1:size(sequence_of_interest,2)
	residuecode = letter2number_map(sequence_of_interest(pos));
	if (residuecode == 0)
		error(['Error: sequence of interest contains undefined residues:' sequence_of_interest]);
	end
	if (residuecode == -1)
		error(['Error: sequence of interest contains ambiguous residues:' sequence_of_interest]);
	end
	if (residuecode > 1)
		focuscolumnlist = [focuscolumnlist pos];
		uniprotoffsetlist = [uniprotoffsetlist next_uniprotoffset];
	end
	if (residuecode == -2 || residuecode > 1)
		next_uniprotoffset = next_uniprotoffset + 1;
	end
end
end

function letter2number_map = create_letter2number_map()
letter2number_map(256) = 0; %initiallize all bytes to 0
letter2number_map('-') = 1;
letter2number_map('A') = 2;
letter2number_map('C') = 3;
letter2number_map('D') = 4;
letter2number_map('E') = 5;
letter2number_map('F') = 6;
letter2number_map('G') = 7;
letter2number_map('H') = 8;
letter2number_map('I') = 9;
letter2number_map('K') = 10;
letter2number_map('L') = 11;
letter2number_map('M') = 12;
letter2number_map('N') = 13;
letter2number_map('P') = 14;
letter2number_map('Q') = 15;
letter2number_map('R') = 16;
letter2number_map('S') = 17;
letter2number_map('T') = 18;
letter2number_map('V') = 19;
letter2number_map('W') = 20;
letter2number_map('Y') = 21;
letter2number_map('B') = -1; %ambiguous : skip sequences containing these
letter2number_map('Z') = -1; %ambiguous : skip sequences containing these
letter2number_map('J') = -1; %ambiguous : skip sequences containing these
letter2number_map('X') = -1; %ambiguous : skip sequences containing these
letter2number_map('U') = -1; %non-standard : skip sequences containing these
letter2number_map('O') = -1; %non-standard : skip sequences containing these
letter2number_map('a') = -2; %non-conserved: skip in seq of interest
letter2number_map('c') = -2; %non-conserved: skip in seq of interest
letter2number_map('d') = -2; %non-conserved: skip in seq of interest
letter2number_map('e') = -2; %non-conserved: skip in seq of interest
letter2number_map('f') = -2; %non-conserved: skip in seq of interest
letter2number_map('g') = -2; %non-conserved: skip in seq of interest
letter2number_map('h') = -2; %non-conserved: skip in seq of interest
letter2number_map('i') = -2; %non-conserved: skip in seq of interest
letter2number_map('k') = -2; %non-conserved: skip in seq of interest
letter2number_map('l') = -2; %non-conserved: skip in seq of interest
letter2number_map('m') = -2; %non-conserved: skip in seq of interest
letter2number_map('n') = -2; %non-conserved: skip in seq of interest
letter2number_map('p') = -2; %non-conserved: skip in seq of interest
letter2number_map('q') = -2; %non-conserved: skip in seq of interest
letter2number_map('r') = -2; %non-conserved: skip in seq of interest
letter2number_map('s') = -2; %non-conserved: skip in seq of interest
letter2number_map('t') = -2; %non-conserved: skip in seq of interest
letter2number_map('v') = -2; %non-conserved: skip in seq of interest
letter2number_map('w') = -2; %non-conserved: skip in seq of interest
letter2number_map('y') = -2; %non-conserved: skip in seq of interest
letter2number_map('b') = -2; %non-conserved: skip in seq of interest
letter2number_map('z') = -2; %non-conserved: skip in seq of interest
letter2number_map('j') = -2; %non-conserved: skip in seq of interest
letter2number_map('x') = -2; %non-conserved: skip in seq of interest
letter2number_map('u') = -2; %non-conserved: skip in seq of interest
letter2number_map('o') = -2; %non-conserved: skip in seq of interest
letter2number_map('.') = -3; %non-conserved: skip in seq of interest, do not advance position
end

function number2letter_map = create_number2letter_map()
number2letter_map(1) = '-';
number2letter_map(2) = 'A';
number2letter_map(3) = 'C';
number2letter_map(4) = 'D';
number2letter_map(5) = 'E';
number2letter_map(6) = 'F';
number2letter_map(7) = 'G';
number2letter_map(8) = 'H';
number2letter_map(9) = 'I';
number2letter_map(10) = 'K';
number2letter_map(11) = 'L';
number2letter_map(12) = 'M';
number2letter_map(13) = 'N';
number2letter_map(14) = 'P';
number2letter_map(15) = 'Q';
number2letter_map(16) = 'R';
number2letter_map(17) = 'S';
number2letter_map(18) = 'T';
number2letter_map(19) = 'V';
number2letter_map(20) = 'W';
number2letter_map(21) = 'Y';
end
