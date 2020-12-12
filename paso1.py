#Step1.py "4D sequence time alignment" 
#Johnatan Quintero V. johnatan2425@hotmail.com

#Step1: Toma archivos A.bvh y B.bvh y genera 3 archivos para cada uno. 
#       Los archivos generados son tipo CSV y contienen la jerarquía,
#       las posiciones globales, y las rotaciones de euler. 

import os                #comandos en bash 

def main():
    commandA1 = "bvh2csv ./A.bvh"      #genera A_rot.csv, A_pos.csv
    commandA2 = "bvh2csv -H ./A.bvh "  #genera A_hierarchy.csv    
    print("Ejecutando el comando: \n", commandA1)
    os.system(commandA1)#escribe en terminal el commandA1
    print("Ejecutando el comando: \n", commandA2)
    os.system(commandA2)#escribe en terminal el commandA2
       
    commandB1 = "bvh2csv ./B.bvh"      #genera B_rot.csv, B_pos.csv
    commandB2 = "bvh2csv -H ./B.bvh "  #genera B_hierarchy.csv
    print("Ejecutando el comando: \n", commandB1)
    os.system(commandB1)#escribe en terminal el commandB1
    print("Ejecutando el comando: \n", commandB2)
    os.system(commandB2)#escribe en terminal el commandB2
main = main()     


