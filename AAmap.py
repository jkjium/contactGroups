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

        self.A2AAA={}
        self.A2AAA['R']='ARG'
        self.A2AAA['H']='HIS'
        self.A2AAA['K']='LYS'
        self.A2AAA['D']='ASP'
        self.A2AAA['E']='GLU'
        self.A2AAA['S']='SER'
        self.A2AAA['T']='THR'
        self.A2AAA['N']='ASN'
        self.A2AAA['Q']='GLN'
        self.A2AAA['C']='CYS'
        self.A2AAA['U']='SEC'
        self.A2AAA['G']='GLY'
        self.A2AAA['P']='PRO'
        self.A2AAA['A']='ALA'
        self.A2AAA['V']='VAL'
        self.A2AAA['I']='ILE'
        self.A2AAA['L']='LEU'
        self.A2AAA['M']='MET'
        self.A2AAA['F']='PHE'
        self.A2AAA['Y']='TYR'
        self.A2AAA['W']='TRP'        
        
        self.A2T={}
        self.A2T['D']='C'
        self.A2T['E']='C'
        self.A2T['H']='C'
        self.A2T['K']='C'
        self.A2T['R']='C'
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