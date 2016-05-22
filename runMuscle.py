import subprocess
import multiprocessing

def runMuscle(muscle_cmd):
    """
    Wywoluje program MUSCLE z wykorzystaniem polecenia wygenerowanego w funkcji generateMuscleCmd(). Korzysta z biblioteki subprocess.
    :param muscle_cmd: polecenie wywolania programu MUSCLE
    :return:
    """
    p = subprocess.Popen(muscle_cmd, shell=True, stderr=subprocess.PIPE)
    p.communicate()

def generateMuscleCmd(klaster, out_name):
    """
    Generuje polecenia wywolania programu MUSCLE
    :param klaster: numer klastra
    :param out_name: nazwa pliku wynikowego
    :return: polecenie wywolania programu MUSCLE
    """
    return "muscle3.8.31_i86win32" + " -in rodziny_homologow\\wynik_" + str(klaster) + ".fasta" + " -fastaout " + out_name

if __name__ == '__main__':

    # Glowna czesc programu, wywolanie MUSCLE
    # program sprawdza liczbe rdzeni i wielowatkowo wywoluje program MUSCLE
    numProcessors = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(numProcessors)

    #Tablica zawierajaca wszystkie wywolania programu MUSCLE
    allMuscleCmd = []

    #Na podstawie danych wynikowych z programu MCL 4207 klastrow
    #Dla kazdego wykonywane jest MSA
    for i in range(1,4208):
        allMuscleCmd.append(generateMuscleCmd(i,"msa_rodzina_" + str(i)))

    i = 0
    tasks = []

    while i < len(allMuscleCmd):
        tasks.append((allMuscleCmd[i],))
        i += 1

    # Run tasks
    results = [pool.apply_async( runMuscle, t ) for t in tasks]

    # Process results
    for i, result in enumerate(results):
        print("Result for cluster %d written to phyi file" % (i))

    pool.close()
    pool.join()
