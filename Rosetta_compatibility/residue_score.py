from pymol import cmd

def bfact_energies(object,pdborig):
        cmd.alter(object, "b=0.0")
        filein = open(pdborig,'r')
        for line in filein:
                if 'pose' in line:
                        break
        k=0
        ener=[]
        for line in filein:
                if 'END_POSE_ENERGIES_TABLE' in line:
                        break
                else:
                        k+=1
                        score = float(line.split()[-5])
                        ener.append(score)
                        cmd.alter('%s and resi %s' %(object,k),'b=%s' %(score))
                        
        filein.close()

        # color the protein based on the new B Factors of the alpha carbons

        cmd.spectrum("b", "blue_white_red", object,minimum=0,maximum=1.4)

cmd.extend("bfact_energies",bfact_energies)

list = cmd.get_object_list('*')
for str in list:
	print(str)
	bfact_energies('%s' %(str),'%s.pdb' %(str)) 
