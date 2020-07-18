'''

The difference between this and the other one will be that for each sample I pick the ethnicity and also two star alleles.
Then at each position I will check the star allele snp file and pull the appropriate genotype from there.
Everything else will be based on the population.

'''




import argparse
import numpy as np
import random
import sys


class GSim(object):
    def __init__(self, freq_file, star_map_file, star_freq_file, as_table=None, count=100, prefix=None, debug=False):
        self.debug = debug
        # Generate a set of variants specified in the manifest to be associated with disease somehow.
        # Specify the frequency file
        self.freq_file = freq_file
        self.star_file = star_map_file
        self.star_freq_file = star_freq_file
        self.as_table_file = as_table

        self.activity_scores = self.read_as_table()

        self.count = int(count)

        self.samples = None

        self.prefix = prefix
        if prefix is None:
            self.prefix = "gsim_%s" % count
        self.run()

    def run(self):

        # Read in the entire star table and exac vcf.  Store in objects
        # star table
        star_map, star_sites = self.make_star_map()

        as_file_samples = self.fetch_as_samples()

        #print(star_map.keys())

        #print(self.activity_scores)

        #exit()


        star_alleles = self.format_stars(as_file_samples)

        # exac
        exac_sites = self.read_exac()

        # get a list of all positions
        all_sites = sorted(set(list(star_sites.keys()) + list(exac_sites.keys())))

        # Read in the sampling frequency for the star alleles


        # create a list of random samples, pass the keys from the star map because some star alleles don't give site info
        self.samples = self.generate_sample_list(star_alleles, star_map.keys(), self.count)

        #print(self.samples)

        if self.as_table_file:
            self.print_activity_scores()


        self.print_header()

        # iterate
        # for s in sites
        for p in all_sites:
            #print(p)
            genotypes = []

            chrom = ""
            pos = ""
            id = ""
            ref = ""
            alts = ""


            # if in stars, just those alleles
            if p in star_sites:

                chrom = star_sites[p]["chrom"]
                pos = star_sites[p]["pos"]
                id = star_sites[p]["id"]
                ref = star_sites[p]["ref"]
                alts = star_sites[p]["alts"]




                for s in self.samples:
                    # check hap 1 and hap 2, generate gt1 and gt2
                    #print(s)
                    hap1 = "CYP2D6_%s_%s" % (s['hap1'], s["hap1_sub"])
                    hap2 = "CYP2D6_%s_%s" % (s['hap2'], s["hap2_sub"])

                    #print(hap1, hap2)

                    if not p in star_map[hap1] or not p in star_map[hap2]:
                        continue
                    gt1 = star_map[hap1][p]
                    gt2 = star_map[hap2][p]
                    genotype = "%s/%s" % (gt1, gt2)
                    genotypes.append(genotype)


            # else, use exac
            elif p in exac_sites:
                #print(exac_sites[p])
                # skip exac indels to avoid overlaps
                indel = False
                if len(exac_sites[p]['ref']) != 1:
                    indel = True
                for a in exac_sites[p]['alts']:
                    if len(a) != 1:
                        indel = True
                #if indel:
                #    continue

                chrom = exac_sites[p]["chrom"]
                pos = exac_sites[p]["pos"]
                id = exac_sites[p]["id"]
                info = exac_sites[p]['info']
                ref = exac_sites[p]['ref']
                alts = exac_sites[p]['alts']
                all_alleles = [ref] + alts

                for s in self.samples:
                    pop = s['pop']

                    allele_1 = self.random_genotype(all_alleles, info, pop)
                    allele_2 = self.random_genotype(all_alleles, info, pop)

                    if allele_1 == ref:
                        gt1 = 0
                    else:
                        gt1 = alts.index(allele_1) + 1
                        # print(gt1)

                    if allele_2 == ref:
                        gt2 = 0
                    else:
                        gt2 = alts.index(allele_2) + 1
                        # print(gt2)

                    genotype = "%s/%s" % (gt1, gt2)
                    genotypes.append(genotype)

            #print(genotypes)

            if len(genotypes) != len(self.samples):
                continue

            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, pos, id, ref, ",".join(alts),
                                                              30, "PASS", ".", "GT", "\t".join(genotypes)))






    def get_stars(self):
        stars = {}
        with open(self.star_freq_file) as f:
            for line in f:
                fields = line.rstrip().split(",")
                stars[fields[0]] = fields[1]

        return stars

    def fetch_as_samples(self):
        samples = []
        with open(self.as_table_file) as f:
            for line in f:
                if line.startswith("Star"):
                    continue
                fields = line.rstrip().split(",")
                samples.append(fields[3])
        return(samples)

    def format_stars(self, stars):

        class SuperStar:
            def __init__(self, star):
                self.star = star
                self.substars = []
                self.activity_score = None

            def random_substar(self):
                return random.choice(self.substars)

            def add_substar(self, substar):
                self.substars.append(substar)

            def print_star(self):
                print("%s, %s, %s" % (self.star, self.substars, self.activity_score))

        star_scores = self.read_as_table()

        star_dict = {}
        for s in stars:
            fields = s.split("_")
            #print(fields)
            superstar = fields[1]
            if superstar not in star_dict:
                star_dict[superstar] = SuperStar(superstar)
                star_dict[superstar].activity_score = star_scores[superstar]
            star_dict[superstar].add_substar(fields[2])

        return star_dict

    def read_as_table(self):
        star_scores = {}
        with open(self.as_table_file) as f:
            for line in f:
                if line.startswith("Star"):
                    continue
                fields = line.rstrip().split(",")
                superstar = fields[1]
                score = fields[6]
                if superstar not in star_scores:
                    star_scores[superstar] = score
        return star_scores

    def make_star_map(self):
        if self.debug:
            print("Making star map", file=sys.stderr)
        star_map = {}
        star_sites = {}
        with open(self.star_file) as f:

            first_line = f.readline()

            #print(first_line.split("\t"))
            #exit()

            stars = first_line.rstrip().split("\t")[9:]
            #print(stars)
            for s in stars:
                star_map[s] = {}
            for line in f:
                fields = line.rstrip().split("\t")


                pos = fields[1]
                star_sites[pos] = {
                    "chrom": fields[0],
                    "pos": fields[1],
                    "id": fields[2],
                    "ref": fields[3],
                    "alts": fields[4].split(",")
                }

                haps = fields[9:]
                for i in range(len(haps)):
                    star_map[stars[i]][pos] = haps[i]
        return(star_map, star_sites)

    def read_exac(self):
        sites = {}
        with open(self.freq_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.rstrip().split()
                chrom = fields[0]
                pos = fields[1]
                id = fields[2]
                ref = fields[3]
                alts = fields[4].split(',')
                info = fields[7]
                sites[pos] = {"chrom": chrom,
                              "pos": pos,
                              "id": id,
                              "ref": ref,
                              "alts": alts,
                              "info": info}
        return(sites)

    def generate_sample_list(self, stars, map_samples, n=100, even=True):
        # generate n samples and assemble them in a list
        #stars = []
        #probs = []
        #for key in star_freqs:
        #    stars.append(key)
        #    probs.append(star_freqs[key])

        samples = []
        if not even:
            for i in range(n):
                samples.append(self.generate_sample(stars, map_samples))
        else:
            for s in [0, 0.5, 1, 1.5, 2]:
                #print("GENERATING! %s" % s)
                for i in range(int(n/5)):
                    samples.append(self.generate_sample(stars, map_samples, s))
        return samples

    def generate_sample(self, stars, map_samples, score=None):
        #Randomly choose an ethnicity and two haplotypes
        sample = {}
        complete = False
        while not complete:
            sample["pop"]= random.choice(["AFR", "AMR", "EAS", "FIN", "NFE"])
            sample["hap1"] = random.choice(list(stars.keys()))
            sample["hap2"] = random.choice(list(stars.keys()))

            sample["hap1_sub"] = stars[sample["hap1"]].random_substar()
            sample["hap2_sub"] = stars[sample["hap2"]].random_substar()

            if score is None:
                return sample
            sample_as = self.get_activity_score(sample["hap1"], sample["hap2"])
            if score == sample_as:
                #print("good score: %s %s" % (sample_as, score))
                return sample
            #print("invalid score! %s %s" % (sample_as, score))



        return sample

    def random_star(self, stars, freqs, map_samples):
        complete = False
        while not complete:
            star = np.random.choice(stars, 1, p=freqs)[0]
            if star in map_samples:
                return star

    def print_header(self):
        samples = []
        count = 0
        for s in self.samples:
            key = "%s_%s_%s" % (s['pop'], s['hap1'].strip("*"), s['hap2'].strip("*"))
            id = self.prefix + "_" + key + "_" + str(count).zfill(5)
            samples.append(id)
            count += 1
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % "\t".join(samples))

    #def get_freqs(self, info):
    #    fields = info.split(";")
    #    allele_freqs = {}
    #    for f in fields:
    #        current = f.split("=")
    #        if current[0] in ["AF_AFR", "AF_AMR", "AF_ASJ", "AF_EAS", "AF_FIN", "AF_NFE", "AF_OTH"]:
    #            freq = current[1].split(',')
    #            total = 0.0
    #            for i in range(0, len(freq)):
    #                if freq[i] == '.':
    #                    freq[i] = 0.0
    #                else:
    #                    freq[i] = float(freq[i])
    #                    total += freq[i]
    #            freq.append(1-total)
    #            allele_freqs[current[0].split("_")[1]] = freq
    #    return allele_freqs

    def find(self, lst, a, b):
        return [i for i, x in enumerate(lst) if x < a or x > b]

    def random_genotype(self, alleles, info, pop):
        # Get alele frequency from info
        freqs = self.get_pop_allele_freq(info, pop)
        # Select allele
        allele = np.random.choice(alleles, p=freqs)
        return allele

    def get_pop_allele_freq(self, info, pop):
        # should be allele count / allele number
        # Split on comma for multi allelic sites
        fields = info.split(";")
        for f in fields:
            if f.startswith("AC_%s" % pop):
                #print("HELLO!")
                subfields = f.split('=')
                allele_counts = subfields[1].split(",")
            if f.startswith("AN_%s" % pop):
                subfields = f.split('=')
                allele_number = int(subfields[1])
        freqs = []
        total = 0
        for c in allele_counts:
            if not allele_number == 0:
                allele_freq = int(c) / allele_number
            else:
                allele_freq = 0
            freqs.append(allele_freq)
            total += allele_freq
        freqs.insert(0, 1-total)
        return freqs

    def get_activity_score(self, hap1, hap2):


        score = float(self.activity_scores[hap1]) + float(self.activity_scores[hap2])
        return score



    def print_activity_scores(self):


        scores = {}

        pops = {}
        file = open('%s.labels.csv' % self.prefix, "w")
        count = 0
        for s in self.samples:
            key = "%s_%s_%s" % (s['pop'], s['hap1'].strip("*"), s['hap2'].strip("*"))
            id = self.prefix + "_" + key + "_" + str(count).zfill(5)
            pop = s['pop']
            hap1 = s["hap1"]
            hap2 = s["hap2"]

            hap1_sub = s["hap1_sub"]
            hap2_sub = s["hap2_sub"]

            score = float(self.activity_scores[hap1]) + float(self.activity_scores[hap2])

            mc = None
            if score == 0:
                mc = "PM"
            elif score == 0.5:
                mc = "IM"
            else:
                mc = "NM"

            file.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (id, pop, hap1, hap1_sub, hap2, hap2_sub, score, mc))
            if score in scores:
                scores[score] += 1
            else:
                scores[score] = 0
            if pop in pops:
                pops[pop] += 1
            else:
                pops[pop] = 0
            count += 1
        #print(scores)
        #print(pops)
        #exit()



"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Generate simulated data based on preset disease patterns')
    parser.add_argument("-f", "--freq_file", help="ExAC frequency file")
    parser.add_argument("-n", "--count", help="Number of samples to generate for each ethnicity")
    parser.add_argument("-p", "--prefix", default=None, help="Output file prefix")
    parser.add_argument("-s", "--star_freqs")
    parser.add_argument("-S", "--star_map")
    parser.add_argument("-a", "--as_table")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()



    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    GSim(options.freq_file, options.star_map, options.star_freqs, options.as_table, options.count, options.prefix, options.debug)

