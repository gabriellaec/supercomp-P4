# -*- coding: utf-8 -*-
import random

i=5
index=0
while i < 200:
     n = i # tamanho da primeira sequência
     m = i # tamanho da segunda sequência
     file = "../relatorio/in_exaustiva_maior/dna{}.seq".format(index) # nome do arquivo a ser gerado
     f = open(file, 'w')
     seq=[str(n)+'\n',
          str(m)+'\n',
          ''.join(random.choices(['A','T','C','G','-'],k=n))+'\n',
          ''.join(random.choices(['A','T','C','G','-'],k=m))]
     f.writelines(seq)
     f.close()
     print(''.join(seq))
     i+=5
     index+=1
