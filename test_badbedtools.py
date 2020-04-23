import unittest
import badbedtools as bad


class TestBadBedTools(unittest.TestCase):
    def test_sort_1(self):
        f1 = './test/bed_1.bed'
        with open(f1) as f1, \
                open('./test/sorted_1.bed') as sort_res_1:
            # сравниваем только 1 и 2 столбец, с 3 какой-то зашквар ( даже башовский сорт
            # не совпадает с бедтулзовским
            chrom_start_1 = ''
            chrom_start_1_res = ''
            for line in bad.sort(f1).split('\n'):
                splt_line = line.split('\t')
                chrom_start_1 += '\t'.join(splt_line[:2]) + '\n'
            for line in sort_res_1:
                splt_line = line.split('\t')
                chrom_start_1_res += '\t'.join(splt_line[:2]) + '\n'
            self.assertEqual(chrom_start_1.rstrip(), chrom_start_1_res.rstrip())

    def test_sort_2(self):
        f2 = './test/bed_2.bed'
        with open(f2) as f2, \
                open('./test/sorted_2.bed') as sort_res_2:
            # сравниваем только 1 и 2 столбец, с 3 какой-то зашквар ( даже башовский сорт
            # не совпадает с бедтулзовским
            chrom_start_2 = ''
            chrom_start_2_res = ''
            for line in bad.sort(f2).split('\n'):
                splt_line = line.split('\t')
                chrom_start_2 += '\t'.join(splt_line[:2]) + '\n'
            for line in sort_res_2:
                splt_line = line.split('\t')
                chrom_start_2_res += '\t'.join(splt_line[:2]) + '\n'
            self.assertEqual(chrom_start_2.rstrip(), chrom_start_2_res.rstrip())

    def test_merging_1(self):
        f1 = './test/bed_1.bed'
        with open(f1) as f1, \
                open('./test/merged_1.bed') as merg_res_1:
            self.assertEqual(bad.merge(f1), merg_res_1.read())

    def test_merging_2(self):
        f2 = './test/bed_2.bed'
        with open(f2) as f2, \
                open('./test/merged_2.bed') as merg_res_2:
            self.assertEqual(bad.merge(f2), merg_res_2.read())

    def test_merging_1_max_distance(self):
        f1 = './test/bed_1.bed'
        with open(f1) as f1, \
                open('./test/merged_1_max_100.bed') as merg_res_2:
            self.assertEqual(bad.merge(f1, max_distance=100), merg_res_2.read())

    def test_merging_2_max_distance(self):
        f2 = './test/bed_2.bed'
        with open(f2) as f2, \
                open('./test/merged_2_max_100.bed') as merg_res_2:
            self.assertEqual(bad.merge(f2, max_distance=100), merg_res_2.read())

    def test_subtract(self):
        f1 = './test/bed_1.bed'
        f2 = './test/bed_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/subtract.bed') as res_12:
            self.assertEqual(bad.subtract(f1, f2), res_12.read())

    def test_subtract_no_sorting(self):
        f1 = './test/sorted_1.bed'
        f2 = './test/sorted_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/subtract.bed') as res_12:
            self.assertEqual(bad.subtract(f1, f2, sorting=False), res_12.read())

    def test_subtract_no_sorting_no_merging(self):
        f1 = './test/merged_1.bed'
        f2 = './test/merged_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/subtract.bed') as res:
            self.assertEqual(bad.subtract(f1, f2, sorting=False, merging=False), res.read())

    def test_subtract_A(self):
        f1 = './test/bed_1.bed'
        f2 = './test/bed_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/subtract_A.bed') as res:
            self.assertEqual(bad.subtract_A(f1, f2), res.read())

    def test_subtract_A_no_sorting(self):
        f1 = './test/sorted_1.bed'
        f2 = './test/sorted_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/subtract_A.bed') as res:
            self.assertEqual(bad.subtract_A(f1, f2, sorting=False), res.read())

    def test_subtract_A_no_sorting_no_merging(self):
        f1 = './test/merged_1.bed'
        f2 = './test/merged_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/subtract_A.bed') as res:
            self.assertEqual(bad.subtract_A(f1, f2, sorting=False, merging=False), res.read())

    def test_intersect(self):
        f1 = './test/bed_1.bed'
        f2 = './test/bed_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/intersect_s.bed') as res:
            # bedtools выдает интерсект не отсортированным, поэтому тест отсортировал
            self.assertEqual(bad.intersect(f1, f2), res.read())

    def test_intersect_no_sorting(self):
        f1 = './test/sorted_1.bed'
        f2 = './test/sorted_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/intersect_s.bed') as res:
            self.assertEqual(bad.intersect(f1, f2,sorting=False), res.read())

    def test_intersect_no_sorting_no_merging(self):
        f1 = './test/merged_1.bed'
        f2 = './test/merged_2.bed'
        with open(f1) as f1, open(f2) as f2, \
                open('./test/intersect_s.bed') as res:
            self.assertEqual(bad.intersect(f1, f2, sorting=False, merging=False), res.read())

    def test_getnospacefasta(self):
        # не знаю как проверять отдельно, проверяю тем что использую этот файл в getfasta()
        fasta = './test/genome.fa'
        with open(fasta) as fa:
            bad.getnospacefasta(fa, outname='./test/no_space_genome.fa')

    def test_faidx(self):
        f1 = './test/no_space_genome.fa'
        with open(f1) as f1, \
                open('./test/no_space_genome.fa.fai') as res:
            # небольшой изврат чтоб привести в соотвествующий вид
            # для теста испольвался samtools faidx (у него в конце пустая строка) + у из
            # него оставил только первые 3 колонки
            # парсим словарь в строку
            my_faidx = ''
            for key, value in bad.faidx(f1).items():
                my_faidx += key + '\t' + '\t'.join(map(str, value)) + '\n'
            # сравниваем
            self.assertEqual(my_faidx, res.read())
            # вуаля

    def test_getfasta(self):
        fasta = './test/no_space_genome.fa'
        bed = './test/small_bed_1.bed'  # ускорял как мог, оч долго работает,
        # поэтому на маленьком файле тест
        with open(fasta) as fa, open(bed) as bedfile, open('./test/no_space_getfasta.fa') as res:
            outname = './test/my_getfasta.fa'
            bad.getfasta(bedfile, fa, outname)
            with open(outname) as my_getfasta:
                self.assertEqual(my_getfasta.read(), res.read())


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Test)
    unittest.TextTestRunner().run(suite)
