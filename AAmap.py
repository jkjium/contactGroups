__all__=['AAmap']

class AAmap(object):# -*- coding: utf-8 -*-

    def __init__(self):
        self.AAA2A={}
        self.AAA2A['ARG']='R'
        self.AAA2A['HIS']='H'
        self.AAA2A['LYS']='K'
        self.AAA2A['ASP']='D'
        self.AAA2A['GLU']='E'
        self.AAA2A['SER']='S'
        self.AAA2A['THR']='T'
        self.AAA2A['ASN']='N'
        self.AAA2A['GLN']='Q'
        self.AAA2A['CYS']='C'
        self.AAA2A['SEC']='U'
        self.AAA2A['GLY']='G'
        self.AAA2A['PRO']='P'
        self.AAA2A['ALA']='A'
        self.AAA2A['VAL']='V'
        self.AAA2A['ILE']='I'
        self.AAA2A['LEU']='L'
        self.AAA2A['MET']='M'
        self.AAA2A['PHE']='F'
        self.AAA2A['TYR']='Y'
        self.AAA2A['TRP']='W'
        
        self.A2T={}
        self.A2T['D']='-'
        self.A2T['E']='-'
        self.A2T['H']='+'
        self.A2T['K']='+'
        self.A2T['R']='+'
        self.A2T['P']='H'
        self.A2T['V']='H'
        self.A2T['M']='H'
        self.A2T['I']='H'
        self.A2T['L']='H'
        self.A2T['F']='H'
        self.A2T['W']='H'
        self.A2T['G']='H'
        self.A2T['A']='H'
        self.A2T['C']='P'
        self.A2T['T']='P'
        self.A2T['Q']='P'
        self.A2T['N']='P'
        self.A2T['Y']='P'
        self.A2T['S']='P'


    def getAAmap(self, AAA):
        return self.AAA2A[AAA]  
    
    def getA2Tmap(self, A):
        return self.A2T[A]