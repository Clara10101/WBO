from Bio import SeqIO
import subprocess, multiprocessing, os, sys
import argparse

def generate_water_cmd(macierz, pliki_fasta_rodzina):
    """
    Generuje polecenia wywolania programu water EMBOSS dla wszystkich sekwencji podanych jako nazwy plikow je zawierajacych
    :param macierz: lokalizacja/nazwa pliku z macierza substytucji PAM/BLOSUM
    :param pliki_fasta_fodzina: lista lokalizacji/nazw plikow z sekwencjami bialkowymi fasta nalezacymi do danej rodziny
    :return: polecenie wywolania programu water
    """
    records = []
    for file in pliki_fasta_rodzina:
        handle = open(file, "rU")
        records.extend(list(SeqIO.parse(handle, "fasta")))
        handle.close()

    from Bio.Emboss.Applications import WaterCommandline
    all_water_cmd = []
    for i in range(len(records)):
        for j in range(len(records)):
            if i < j:

                water_cmd = WaterCommandline(gapopen=100, gapextend=10)#maksymalne wartosci aby uzyskac uliniowienia bezspacjowe
                water_cmd.asequence = "asis:" + str(records[i].seq)
                water_cmd.bsequence = "asis:" + str(records[j].seq)
                water_cmd.stdout = True
                water_cmd.sprotein=True
                water_cmd.datafile=macierz
                all_water_cmd.append(str(water_cmd))

    return all_water_cmd

def parseOutput(waterOutput):
    '''
    Parsuje wynik programu water, wybierajac 28 linijke, w ktorej podana jest wartosc score
    :param waterOutput: wynik programu water
    :return: score alignmentu logalnego jako float
    '''
    return float(waterOutput.split('\r\n')[28].split()[2])

def runWater(water_cmd):
    '''
    Wywoluje program water EMBOSS z wykorzystaniem polecenia wygenerowanego w funkcji generate_water_cmd(). Korzysta z biblioteki subprocess
    :param water_cmd: polecenie wywolania programu water
    :return: score alignmentu lokalnego
    '''
    p = subprocess.Popen(water_cmd, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)#, stderr=subprocess.PIPE)
    alignment_result = p.communicate(os.linesep)[0].rstrip()
    return parseOutput(alignment_result)

def runAsynchronously(macierz,rodzina):
    '''
    Wywoluje program water asynchronicznie dla wszystkich sekwencji w rodzinie i dla podanej w parametrze macierzy
    :param macierz: lokalizacja/nazwa pliku z macierza substytucji PAM/BLOSUM
    :param rodzina: lista lokalizacji/nazw plikow z sekwencjami bialkowymi fasta nalezacymi do danej rodziny
    :return: sredni score dla sekwencji podanych w liscie
    '''

    # Glowna czesc programu, wywolanie water z wykorzystaniem biblioteki multiprocessing,
    # program sprawdza liczbe rdzeni i wielowatkowo wywoluje program water
    numProcessors = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numProcessors)

    all_water_cmd = generate_water_cmd(macierz, rodzina)

    i = 0
    tasks = []

    while i < len(all_water_cmd):
        tasks.append((all_water_cmd[i],))
        i += 1

    # Run tasks
    results = [pool.apply_async(runWater, t) for t in tasks]

    scoreSum = 0
    # Process results
    for i, result in enumerate(results):
        score = result.get()
        scoreSum += score

    pool.close()
    pool.join()

    return scoreSum / len(results)

if __name__ == '__main__':

    #glowna czesc programu
    #parametry podawane w linii polecen
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", help="nazwa katalogu zawierajacego macierze substytucji PAM/BLOSUM", type=str, required=True)
    parser.add_argument("-f", "--fasta", help="katalog zawierajacy sekwencje w formacie fasta", type=str)
    parser.add_argument('rest', nargs=argparse.REMAINDER)#pozostale argumenty to nazwy plikow w formacie fasta
    args = parser.parse_args()

    matrix_directory = args.directory

    matrixFiles = []

    for file in os.listdir(matrix_directory):
        #wykorzystywane macierze substytucji sa zapisane w formacie txt
        if file.endswith(".txt"):
            matrixFiles.append(file)

    #katalog z sekwencjami w formacie fasta to parametr opcjonalny
    #alternatywnie moga zostac podane bezposrednio nazwy plikow fasta
    if args.fasta:
        family_directory = args.fasta

        familyFiles = []

        #zakladamy ze sekwencje nalezace do rodziny homologow zapisane sa w katalogu o podanej nazwie
        #sa to wszystkie pliki fasta w pliku

        for file in os.listdir(family_directory):
            #zalozenie ze w folderze sa jedynie pliki zapisane w formacie fasta
            #w przypadku innego formatu nie zostana one sparsowane przez wykorzystywany modul biopythona
            familyFiles.append(os.path.join(family_directory,file))

    else:

        # wyjatek jesli nie zostal podany katalog z sekwencjami ani pliki fasta do analizy
        try:
            args.rest
        except AttributeError:
            print "Nalezy podac katalog z sekwencjami lub nazwy plikow fasta do analizy"
            sys.exit(1)

        familyFiles = args.rest

    scores = []

    #wywolanie dla kazdej macierzy w folderze i zliczenie wynikow
    for matrix in matrixFiles:
        score = runAsynchronously(os.path.join(matrix_directory,matrix), familyFiles)
        scores.append(score)
        print '\nSredni koszt uliniowienia dla macierzy ' + matrix[:-4] + ' wynosi ' + str(score)

    print 'Macierz generujaca najlepsze bezspacjowe uliniowienie lokalne dla danej rodziny homologow to macierz ' + matrixFiles[scores.index(max(scores))][:-4]
    print 'Wynik uliniowienia wynosi ' + str(max(scores))
