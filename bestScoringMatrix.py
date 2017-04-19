import os
import subprocess
import datetime


def get_immediate_subdirectories(a_dir):
    return [os.path.join(a_dir, name) for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def bestScoringMatrix(output_file):
    '''
    Funkcja pomocnicza, niewykorzystywana w glownym programie do wyznaczania najlepszej macerzy substytucji
    PAM/BLOSUM dla kazdej z rodzin sekwencji homologicznych
    :return:
    '''

    resultFile = open(output_file,'w')

    #zmienne z folderami zawierajacymi wszystkie rodziny
    #sekwencje fasta nalezace do rodziny x archea znajduja sie w katalogu archea/rodzina_x
    #analogicznie dla bakterii, wirusow i drozdzy

    archea_families = get_immediate_subdirectories('archea')
    bacteria_families = get_immediate_subdirectories('bakterie')
    yeast_families = get_immediate_subdirectories('drozdze')
    viruses_families = get_immediate_subdirectories('wirusy')

    len_archea_families = len(archea_families)
    len_bacteria_families = len(bacteria_families)
    len_yeast_families = len(yeast_families)
    len_viruses_families = len(viruses_families)

    all_families_directories = archea_families + bacteria_families + yeast_families + viruses_families
    #print all_families_directories potrzebne !!
    #macierze PAM/BLOSUM zapisane odpowiednio w folderach pam i blosum

    print datetime.datetime.now().time()

    #dla kazdej rodziny w grupie wyznaczamy najlepsze macierze PAM i BLOSUM
    for i,directory in enumerate(all_families_directories):

        #dla skrocenia czasu obliczen - co 10 rodzina
        if i % 10 == 0:

            print 'Wyznaczanie najlepszej macierzy pam'

            cmd = 'python skrypt.py -d pam -f ' + directory
            #pobranie wyniku wywolania
            result = subprocess.check_output(cmd, shell=True).split()

            bestPAM = result[-5]
            bestPAMScore = result[-1]

            print 'Wyznaczenie najlepszej macierzy blosum'

            cmd = 'python skrypt.py -d blosum -f ' + directory
            result = subprocess.check_output(cmd, shell=True).split()

            bestBLOSUM = result[-5]
            bestBLOSUMScore = result[-1]

            #zapis danych wynikowych w formacie csv
            if i < len_archea_families:
                resultFile.write('archea\t')
            elif i < len_archea_families + len_bacteria_families:
                resultFile.write('bacterie\t')
            elif i < len_archea_families + len_bacteria_families + len_yeast_families:
                resultFile.write('drozdze\t')
            else:
                resultFile.write('wirusy\t')

            resultFile.write(directory.split('\\')[1]+'\t'+ bestPAM + '\t' + bestPAMScore + '\t'+ bestBLOSUM + '\t' + bestBLOSUMScore +'\n')

    print datetime.datetime.now().time()

bestScoringMatrix('wyniki.csv')
