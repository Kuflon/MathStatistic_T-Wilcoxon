import sys
from copy import copy

before = [25, 50, 77, 76, 67, 75, 77, 71, 63, 98, 60]
after =  [64, 77, 74, 95, 105, 83, 73, 75, 101, 97, 78]


swapper = lambda x, y: (copy(y), copy(x))


def check_n(n1, n2):
    '''
    Checks for the number of items in each selection.
    -------------------------------------------------
    Проверяет количество элементов в каждой выборке.
    '''
    if(n1 != n2):
        sys.exit("The size of the first sample is not equal to the size of the second!")
    if(n1 < 5):
        sys.exit("The samples are too small!")
    if(n1 > 50):
        sys.exit("Too large a sample n must be < 51 !")


def t_table(n):
    tbl005 =   [0, 2, 3, 5, 8, 10, 13, 17, 21, 25, 30, 35, 41, 47, 53, 60, 67, 75, 83, 91, 100, 110, 119, 130, 140, 151, 163, 175, 187, 200, 213, 227, 241, 256, 271, 286, 302, 319, 336, 353, 371, 389, 407, 426, 446, 466]

    tbl001 =    [float('nan'), float('nan'), 0, 1, 3, 5, 7, 9, 12, 15, 19, 23, 27, 32, 37, 43, 49, 55, 62, 69, 76, 84, 92, 101, 110, 120, 130, 140, 151, 162, 173, 185, 198, 211, 224, 238, 252, 266, 281, 296, 312, 328, 345, 362, 379, 397]

    g005 = tbl005[n - 5]
    g001 = tbl001[n - 5]

    return g005, g001


def t_wilcoxon(bef, aft):
    '''
    Hypotheses:
        H0: the Intensity of changes in a characteristic direction does not exceed the intensity changes in an atypical direction
        H1: the Intensity shifts in the direction of typical superior to intensity changes in an atypical direction
    -----------------------------------------
    Гипотезы:
        H0: Интенсивность сжвтгов в типичном направлении не превосходит интенсивности сдвигов в нетипичном направлении
        H1: Интенсивность сжвтгов в типичном направлении превосходит интенсивности сдвигов в нетипичном направлении
    '''
    check_n(len(bef), len(aft))

    shift = []
    absshift = []
    temp = 0
    for i in range(len(bef)):                       # Step 1.
        shift.append(aft[i - temp] - bef[i - temp]) # Сalculate the shifts
        absshift.append(abs(shift[i]))

        if shift[i - temp] == 0:                    # Step 2.
            shift.pop(i - temp)                     # Exclude zero shifts
            aft.pop(i - temp)                       # ------------------------
            bef.pop(i - temp)                       # Исключить нулевые сдвиги
            absshift.pop(i - temp)                  #
            temp += 1                               #

    if(len(shift) < 5):
        sys.exit("After excluding zero reactions, the number of observations does not meet the restrictions. There is no way to continue the test.\nПосле исключения нулевых реакций, количество наблюдений не удовлетворяет ограничениям. Нет возможности продолжить тест.")

    t = plus = minus = 0
    difference = []
    for i in range(len(shift)):         # Step 3.
        if shift[i] < abs(shift[i]):    # Determine the prevailing shifts(typical)
            difference.append('-')      #
            minus += 1                  # -----------------------------------------
        else:                           # Определить преобладающие сдвиги(типичные)
            difference.append('+')      #
            plus += 1                   #

    if plus > minus:
        res = '-'
    else:
        res = '+'

    for j in range(len(shift)):                                                                # Step 4.
        for i in range(len(shift) - 1):                                                        # Sort all values in ascending order
            if(absshift[i] > absshift[i + 1]):                                                 # ------------------------------------------
                absshift[i], absshift[i + 1] = swapper(absshift[i], absshift[i + 1])           # Отсортировать все значения по возрастанию
                difference[i], difference[i + 1] = swapper(difference[i], difference[i + 1])   #

    ranks = []
    for i in range(len(absshift)):    # Step 5. Assign ranks to values
        ranks.append(i + 1)           # Присвоить значениям ранги

    absshift.append(0) #blank
    same_ranks = []
    for i in range(len(absshift) - 2):        # Step 6.
        while(absshift[i] == absshift[i+1]):  # Make an adjustment for repeated values
            same_ranks.append(ranks[i])       #
            i += 1                            # ------------------------------------------
        if(absshift[i] == absshift[i-1]):     #
            same_ranks.append(ranks[i])       # Сделать поправку на повторяющиеся значения


        if same_ranks:
            sum = 0
            for k in same_ranks:
                sum += k
            for j in range(len(same_ranks)-1, -1, -1): # Assign a rank equal to the arithmetic mean of the corresponding repeating ranks
                ranks[i - j] = sum / len(same_ranks)   # Присвоить ранг, равный среднему арифметическому соответствующих повторяющихся рангов
            same_ranks = []

    for i in range(len(ranks)):
        if res == '+' and difference[i] == '+':
            t += ranks[i]
        elif res == '-' and difference[i] == '-':
            t += ranks[i]

    t005, t001 = t_table(len(shift))

    print(f'\nt = {t}, t_005 = {t005}, t_001 = {t001}\n')
    if(t > t005):
        print('t > t_005 there is no reason to reject Ho with a confidence level of 95%. The shift in the typical direction is not significant\n')
        if(t > t001):
            print('t > t_001 there is no reason to reject Ho with a confidence level of 99%. The shift in the typical direction is not significant\n')
    else:
        print('t < t_table. H1. The shift in the typical direction is significant \n')

t_wilcoxon(before, after)
