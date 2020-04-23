# Dependencies
from operator import itemgetter
import locale
# check file extension
def check_extension(extension: str):
    bed_variant_extensions = ['bed' + str(num) for num in range(1, 13)] + ['bed']
    extension_variations = ['vcf', 'gff'] + bed_variant_extensions
    if extension not in extension_variations:
        raise TypeError('Need file in bed/vcf/gff format for input')  # можно ввести проверку формата


def split_comments(bed_vcf_gff, extension):
    comment_lines_list = []
    splt_lines_list = []
    for line in bed_vcf_gff:
        if line.startswith('#'):
            comment_lines_list.append(line.rstrip())
        else:
            splt_line = line.rstrip().split('\t')
            if extension in ['bed' + str(num) for num in range(1, 13)] + ['bed']:
                splt_line[1] = int(splt_line[1])
                splt_line[2] = int(splt_line[2])
            elif extension == 'vcf':
                splt_line[1] = int(splt_line[1])
            elif extension == 'gff':
                splt_line[3] = int(splt_line[3])
                splt_line[4] = int(splt_line[4])
            splt_lines_list.append(splt_line)

    return [comment_lines_list, splt_lines_list]


def extract_intervals(intervals_list: list, extension):
    just_intervals = []

    if extension in ['bed' + str(num) for num in range(1, 13)] + ['bed']:
        for line in intervals_list:
                just_intervals += [[line[0], line[1], line[2]]]
    elif extension == 'vcf':
        for line in intervals_list:
                just_intervals += [[line[0], line[1], line[1] + 1]]
    elif extension == 'gff':
        for line in intervals_list:
                just_intervals += [[line[0], line[3], line[4]]] #??? в gff точно так?

    return just_intervals


def overlap(interval_1, interval_2, max_distance=0):
    chrom1, start1, end1 = interval_1
    chrom2, start2, end2 = interval_2
    # [s,e)  а в gff и в vcf наверно входит? Проверить
    if chrom1 == chrom2:
        if end1 + max_distance >= start2:  # разбиваем условия для ускорения программы
            if start1 - max_distance <= end2:
                return True
        # ТУТ для мерджа сделал >= <=, бедтулс объединяет состыкующиеся интервалы
        # но технически они не перекрываются!!! Так что название фукции не очень, но оставлю пока

    return False


def delete_overlaps_intervals_in_list(intervals_list, intervals_list_sub,
             extension, extension_sub):
    # make undesirable_intervals
    intervals = extract_intervals(intervals_list, extension)
    intervals_sub = extract_intervals(intervals_list_sub, extension_sub)

    no_overlap_intervals = []
    count_overlaps = 0
    for interval in intervals:
        overlaping = False
        for subtracted_interval in intervals_sub:
            if interval[0] != subtracted_interval[0]:
                continue
            if interval[2] < subtracted_interval[1]:
                continue
            if interval[1] > subtracted_interval[2]:
                continue

            overlaping = True  # skip overlapped interval
            count_overlaps += 1
        if not overlaping:
            no_overlap_intervals.append(interval)

    return [no_overlap_intervals,
            count_overlaps,
            len(no_overlap_intervals)]


def merge_interval(interval_1, interval_2):
    if interval_1[0] != interval_2[0]:
        raise ValueError('Not identical chromosome names for merge two interval')
    chrom = interval_1[0]
    start_1, end_1 = map(int, interval_1[1:])
    start_2, end_2 = map(int, interval_2[1:])
    start = min(start_1, start_2)
    end = max(end_1, end_2)
    return [chrom, start, end]


def merge_intervals_list(intervals_list: list, extension, max_distance=0):
    # extract intervals for merge
    intervals = extract_intervals(intervals_list, extension)
    merged_intervals = []
    count_merges = 0

    for i, interval in enumerate(intervals):
        if i == 0:
            previous_interval = interval
            continue
        if overlap(interval, previous_interval, max_distance):
            previous_interval = merge_interval(interval, previous_interval)
            count_merges += 1
            if i == len(intervals) - 1:  # Если последний интервал перекрывается с предыдущим записываем его
                merged_intervals.append(previous_interval)
        elif i != len(intervals) - 1:
            merged_intervals.append(previous_interval)
            previous_interval = interval
        else:
            merged_intervals.append(previous_interval)
            merged_intervals.append(interval)

    return [merged_intervals, count_merges]


def intersect_intervals(intervals_1, intervals_2, subtract=False):
    chromosomes_2 = set(map(lambda x: x[0], intervals_1))
    intersected_intervals = []
    for interval_1 in intervals_1:

        chrom1, start1, end1 = interval_1

        if chrom1 not in chromosomes_2:
            if subtract:
                intersected_intervals.append([chrom1, start1, end1])
            continue

        full_overlaped = False
        for interval_2 in intervals_2:
            chrom2, start2, end2 = interval_2
            # скипаем все ненужное

            if (chrom1 == chrom2) and (end1 > start2) and (start1 < end2):
                # разбираем все случаи
                #  ----
                # ------
                if (start1 >= start2) and (end1 <= end2):
                    # тут нащ интервал вычтен, просто прерываем цикл
                    if subtract:
                        full_overlaped = True
                        break
                    else:
                        intersected_intervals.append([chrom1, start1, end1])
                        break  # сохраняем исходные start1 и end1
                #   -----
                # -----
                elif (start1 >= start2) and (start1 < end2):
                    if subtract:
                        start1 = end2
                        end1 = end1
                    else:
                        intersected_intervals.append([chrom1, start1, end2])
                        start1 = end2  # end1 сохраняется
                # -----
                #   -----
                elif (end1 > start2) and (end1 <= end2):
                    if subtract:
                        intersected_intervals.append([chrom1, start1, start2])
                        full_overlaped = True
                        break
                    else:
                        intersected_intervals.append([chrom1, start2, end1])
                        break
                # ------
                #  ----
                elif (start2 > start1) and (end2 < end1):
                    if subtract:
                        # записываем левый интервал, он уже ни с чем не пересечется
                        intersected_intervals.append([chrom1, start1, start2])
                        # с правым работаем дальше
                        start1 = end2
                        end1 = end1
                    else:
                        intersected_intervals.append([chrom1, start2, end2])
                        start1 = end2  # end1 остается


            # добавляем интервал новый если чет осталось
        if not full_overlaped and subtract:  # Если не перекрылся ни с кем при вычитании добавляем
            intersected_intervals.append([chrom1, start1, end1])
    return intersected_intervals


def intersect_intervals_list(intervals_list, intervals_list_int, extension, extension_int,
              subtract=False):
    # Если использовать как сейчас лишний шаг, но если опционально выключать сорт
    # и мердж лучше оставить
    intervals = extract_intervals(intervals_list, extension)
    intervals_int = extract_intervals(intervals_list_int, extension_int)

    result_intervals = intersect_intervals(intervals, intervals_int, subtract)

    return result_intervals


# starts from 0 possitions
# end  non-inclusive


def sort_intervals_list(intervals_list: list, extension):  # return sorted str of lines from file
    if extension in ['bed' + str(num) for num in range(1, 13)] + ['bed']:
        intervals_list.sort(key=itemgetter(0, 1))
    elif extension == 'vcf':
        intervals_list.sort(key=itemgetter(0, 1))
    elif extension == 'gff':
        intervals_list.sort(key=itemgetter(0, 3))

    return intervals_list


def convert_intervals_to_str(intervals, comments=[]):
    result = '\n'.join(comments)
    for line in intervals:
        # Сделаем все стркоми
        strline = map(str, line)
        joinline = '\t'.join(strline)
        result += joinline + '\n'
    return result



# Sort bed/vcf/gff
def sort(bed_vcf_gff_file, extension='bed'):
    """
    Order the intervals in a file.
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format
    :param extension: <bed/gff/vcf>
    :return: string
    """
    check_extension(extension)
    # Сплитим на комменты и интервалы + по колонкам
    comment_lines_list, intervals_list = split_comments(bed_vcf_gff_file, extension)
    # Сортируем интервальный лист по колонкам ссоотвестующим расширению
    sorted_intervals_list = sort_intervals_list(intervals_list, extension)
    # Соединяем комменты и интервалы отсортированные
    results = convert_intervals_to_str(sorted_intervals_list, comment_lines_list)
    # возвращаем строку готовую для записи в файл
    return results

def subtract_A(bed_vcf_gff_file,bed_vcf_gff_subtracted_file,
             extension_1='bed', extension_2='bed', sorting=True, merging=True):
    """
        Delete intervals that is overlapped by another feature(s)
        :param bed_vcf_gff_file: File in <bed/gff/vcf> format.
        :param bed_vcf_gff_subtracted_file: Subracted file in <bed/gff/vcf> format.
        :param extension_1: <bed/gff/vcf>
        :param extension_2: <bed/gff/vcf>
        :param sorting: Add sorting step. Default = True
        :param merging: Add merging step. Default = True
        :return: string
        """
    # чекаем расширение
    check_extension(extension_1)
    check_extension(extension_2)
    # Сплитим
    comment_lines_list, intervals_list = split_comments(bed_vcf_gff_file,
                                                        extension_1)
    comment_lines_list_sub, intervals_list_sub = split_comments(bed_vcf_gff_subtracted_file,
                                                                extension_2)
    if sorting:
        # Сортируем
        intervals_list = sort_intervals_list(intervals_list, extension_1)
        intervals_list_sub = sort_intervals_list(intervals_list_sub, extension_2)
    if merging:
        # Мерджим
        intervals_list, count_merges = merge_intervals_list(intervals_list, extension_1, )
        intervals_list_sub, count_merges_sub = merge_intervals_list(intervals_list_sub,
                                                                    extension_2)
        # После merge у нас уже все в bed формате
        extension_1 = 'bed'
        extension_2 = 'bed'

    # Вычитаем
    no_overlap_interval, count_overlaps, count_no_overlaps = delete_overlaps_intervals_in_list(intervals_list,
                                                                                      intervals_list_sub,
                                                                                      extension_1,
                                                                                      extension_2)
    # Можно добавить добавление комментария о вычитании файла
    result = convert_intervals_to_str(no_overlap_interval, comment_lines_list)
    # Статистику выводим (для консольной версии)
    # print(f"""There are {count_overlaps} overlaps in file1 with file2.
    # And {count_no_overlaps} are not overlaped intervals.""")
    return result
# Subtract intrevals bed/vcf/gff from bed/vcf/gff
def subtract(bed_vcf_gff_file,bed_vcf_gff_subtracted_file,
             extension_1='bed', extension_2='bed', sorting=True, merging=True):
    """
    Removes the portion(s) of an interval that is overlapped by another feature(s)
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format.
    :param bed_vcf_gff_subtracted_file: Subracted file in <bed/gff/vcf> format.
    :param extension_1: <bed/gff/vcf>
    :param extension_2: <bed/gff/vcf>
    :param sorting: Add sorting step. Default = True
    :param merging: Add merging step. Default = True
    :return: string
    """
    subtract = True
    # чекаем расширение
    check_extension(extension_1)
    check_extension(extension_2)
    # Сплитим
    comment_lines_list, intervals_list = split_comments(bed_vcf_gff_file,
                                                        extension_1)
    comment_lines_list_sub, intervals_list_sub = split_comments(bed_vcf_gff_subtracted_file,
                                                                extension_2)
    if sorting:
        # Сортируем
        intervals_list = sort_intervals_list(intervals_list, extension_1)
        intervals_list_sub = sort_intervals_list(intervals_list_sub, extension_2)

    if merging:
        # Мерджим
        intervals_list, count_merges = merge_intervals_list(intervals_list, extension_1, )
        intervals_list_sub, count_merges_sub = merge_intervals_list(intervals_list_sub,
                                                   extension_2)
        # После merge у нас уже все в bed формате
        extension_1 = 'bed'
        extension_2 = 'bed'



    # Вычитаем
    no_overlap_interval = intersect_intervals_list(intervals_list,
                                                  intervals_list_sub,
                                                  extension_1,
                                                  extension_2,
                                                  subtract=True)
    # Можно добавить добавление комментария о вычитании файла
    result = convert_intervals_to_str(no_overlap_interval, comment_lines_list)
    # Статистику выводим (для консольной версии)
    # print(f"""There are {count_overlaps} overlaps in file1 with file2.
    # And {count_no_overlaps} are not overlaped intervals.""")
    return result




def merge(bed_vcf_gff_file, extension='bed', max_distance=0, sorting=True):
    """
    Merges overlapping sorted BED/GFF/VCF entries into a single interval.
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format.
    :param extension: <bed/gff/vcf>.
    :param max_distance: Maximum distance between features allowed for features to be merged. Default is 0
    :param sorting: Add sorting step. default = True
    :return: string
    """
    # Сплитим
    comment_lines_list, intervals_list = split_comments(bed_vcf_gff_file, extension)
    if sorting:
        # Сортируем
        intervals_list = sort_intervals_list(intervals_list, extension)
    # Мерджим
    merged_intervals, count_merges = merge_intervals_list(intervals_list,
                                                          extension,
                                                          max_distance)
    # тут на выходе получаем bed файл всегда, поэтому комменты обнуляем
    # но потом можно будет сделать мердж с выходом в других форматах
    comment_lines_list = []
    result = convert_intervals_to_str(merged_intervals, comment_lines_list)
    # Возвращаем строку ! (тут проблема с инпутом аутпутом для разных функций,
    # надо каждый раз строку превращать в файл io.StrigIO(),
    # Можно потом добавить вараинт инпута в виде строки, будет удобнее
    # Статистика (для консольной версии)
    # print(f"""Merged overlaps with max distance: {max_distance}.
    # There are {count_merges} overlaps merged.""")
    return result

# intersect == subtract
# reverse intersect == intersect
def intersect(bed_vcf_gff_file, bed_vcf_gff_intersect_file,
              extension='bed', extension_int='bed', subtract=False,
              sorting=True, merging=True):
    """
    Report overlaps between two feature files.
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format
    :param bed_vcf_gff_intersect_file: File in <bed/gff/vcf> format
    :param extension: <bed/gff/vcf>
    :param extension_int: <bed/gff/vcf>
    :param subtract: Return non intersected regions
    :param sorting: Add sorting step. default = True
    :param merging: Add merging step. default = True
    :return: string
    """
    # Проверяем расширения
    check_extension(extension)
    check_extension(extension_int)

    # Сплитим
    comment_lines_list, intervals_list = split_comments(bed_vcf_gff_file, extension)
    comment_lines_list_int, intervals_list_int = split_comments(bed_vcf_gff_intersect_file,
                                                                extension_int)
    if sorting:
        # Сортируем можно сделать ОПЦИОНАЛЬНО для ускорения
        intervals_list = sort_intervals_list(intervals_list, extension)
        intervals_list_int = sort_intervals_list(intervals_list_int, extension_int)
    if merging:
        # Мерджим
        intervals_list, count_merges = merge_intervals_list(intervals_list,
                                                            extension)
        intervals_list_int, count_merges_int = merge_intervals_list(intervals_list_int,
                                                       extension_int)
        # После merge у нас уже все в bed формате
        extension = 'bed'
        extension_int = 'bed'

    # Интерсектим интервалы
    intersected_intervals = intersect_intervals_list(intervals_list,
                                                     intervals_list_int,
                                                     extension,
                                                     extension_int,
                                                     subtract)
    # Конвертим листы интервалов в строку
    result = convert_intervals_to_str(intersected_intervals)
    return result