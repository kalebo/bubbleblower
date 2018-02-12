#!/usr/bin/env python3
import random
from enum import Enum
from math import e
import ipdb

"""
TODO:
  [ ]  add merge capability
  [ ]  add graph generation
  [ ]  
"""

class KmerTracker():
    def __init__(self, k):
        self.used_kmers = set()
        self.k = k

    def isUsed(self, kmer) -> bool:
        return kmer in self.used_kmers

    def add(self, kmer):
        assert len(kmer) == self.k
        self.used_kmers.add(kmer)

class KmerGenerator():
    alphabet = "ACGT"
    max_tries = len(alphabet)

    def __init__(self, tracker, kmer_len=9):
        self._bootstrap(kmer_len)
        self.tracker = tracker

    def _bootstrap(self, k):
        self.seq = [random.choice(self.alphabet) for i in range(k)]


    def _find_unique_kmer(self, seq):
        """Find a unique version of seq (which should be len = k-1) with a random base"""
        for _ in range(self.max_tries):
            candidate = random.choice(self.alphabet)
            candidate_kmer = ''.join(seq) + candidate
            if self.tracker.isUsed(candidate_kmer) is False:
                return candidate, candidate_kmer

        assert False, "Ran out of unique kmers!"

    def random_mutation(self, base_seq: str):
        seq = base_seq[1:]
        last_base, kmer = self._find_unique_kmer(seq)
        self.tracker.add(kmer)
        return last_base

    def generate(self) -> str:
        while True:

            self.seq.pop(0)
            last_base, kmer = self._find_unique_kmer(self.seq)
            self.seq.append(last_base)

            self.tracker.add(kmer)

            yield kmer

class KmerManager():
    def __init__(self, k):
        self.k = k
        self.tracker = KmerTracker(k)
        self.generator = KmerGenerator(self.tracker, k)

class Mode(Enum):
    UNSEEDED = 0 # Needs initial kmer seed
    READY = 1    # Ready to merge or diverge
    COOLDOWN = 2 # Recently diverged; needs to cool down before READY
    PAUSED = 3   # Not extending kmers at all to create indel blocks

class SpeciesContainer():

    def __init__(self, species_id, k): # Cooldown time is k extensions
        self.mode = Mode.UNSEEDED
        self.seq = ""
        self.sid = species_id
        self.k = k
        self.cooldown_timer = 0
        self.last_mutation_timer = 0

    def extend(self, seq, diverge=False): # potentially a simple, divergent, or convergent extend
        if self.mode == Mode.UNSEEDED:
            self.seq = seq
            self.mode = Mode.READY

        elif self.mode == Mode.PAUSED:
            assert False, "Not implemented."

        else: # READY or COOLDOWN
            self.seq += seq[-1]
            if diverge:
                self.cooldown_timer = self.k
                self.last_mutation_timer = 0
                self.mode = Mode.COOLDOWN
            elif self.mode == Mode.COOLDOWN:
                self.cooldown_timer -= 1
                if self.cooldown_timer == 0:
                    self.mode = Mode.READY

            self.last_mutation_timer += 1

    def tip_kmer(self) -> str:
        """Returns the kmer at the tip of the sequence"""
        return self.seq[-self.k:]

class SpeciesManager():
    def __init__(self, species_count, kmer_manager, allow_cross_merges=False):
        self.manager = kmer_manager
        self.kmer_generator = kmer_manager.generator
        self.generator = kmer_manager.generator.generate()
        self.tracker = kmer_manager.tracker
        self.cross_merges = allow_cross_merges
        self.iteration = 0
        self.mu = 1/16 # mutation/iteration
        k = kmer_manager.k

        self.species = { i : SpeciesContainer(i, k) for i in range(species_count) }

        # Initially all species will be present in the first kmer and thus a single group
        self.groups = [ frozenset(i for i in range(species_count)) ]
        self.prev_groups = {}

    def _split_group(self, group) -> (frozenset, frozenset):
        if len(group) <= 1: # This avoids bulges
            assert False, "Will not split groups of 1 as this creates \"bulges\""

        right = set()
        left = set()
        for i in group:
            if random.random() >.5:
                right.add(i)
            else:
                left.add(i)

        # make sure that there is at least one species in each group
        if len(right) == 0:
            right.add(left.pop())
        elif len(left) == 0:
            left.add(right.pop())
        assert len(left) !=0 and len(right) != 0

        fright = frozenset(right)
        fleft = frozenset(left)

        self.prev_groups[fright] = group
        self.prev_groups[fleft] = group

        self.groups.remove(group)
        self.groups.append(fright)
        self.groups.append(fleft)

        return fleft, fright

    def _is_group_mergeable(self, union):
        for i in union:
            if self.species[i].mode != Mode.READY:
                return False
        return True

    def _merge_groups(self, left, right) -> (bool, frozenset):
        group = left.union(right)
        if group == self.prev_groups[left] and \
           group == self.prev_groups[right] and \
           self._is_group_mergeable(group):

            self.prev_groups.pop(right)
            self.prev_groups.pop(left)

            self.groups.remove(left)
            self.groups.remove(right)
            self.groups.append(group)

            return True, group

        return False, {}

    def _p_mutation(self, species: SpeciesContainer):
        """
        Models a poisson process for the probability of a divergent mutation.
        1 - e^(-µx)
            µ is average rate of mutation for a species which lessens as iterations increase past a threshold.
            x is how many iterations have passed since a mutation in a species.
        """
        return 1 - e**(-(self.mu if self.iteration < 25 else 0) * species.last_mutation_timer) #todo: scale mu instead of nuking it.

    def _find_containing_group(self, sid):
        for g in self.groups:
            if sid in g:
                return g

    def iterate(self):
        # try merges
        # Diverging groups are defined by having a different base extensions, all other tips get the non-divergent base
        # A merge requires both groups to have their union equal both of their previous group sets.
        kmer = next(self.generator)

        diverging_species = set()
        for species in self.species.values():
            if random.random() < self._p_mutation(species):
                diverging_species.add(species.sid)

        diverging_species = list(diverging_species) # need to get a data structure that can be index to randomized
        diverging_species_groups = set()
        random.shuffle(diverging_species)
        for sid in diverging_species:
            group = self._find_containing_group(sid)
            if len(group) > 1:
                # split and then find which other species would have diverged together.
                left, right = self._split_group(group)
                if sid in left:
                    diverging_species_groups.add(left)
                else:
                    diverging_species_groups.add(right)

        for group in self.groups:
            if group in diverging_species_groups:
                # Generate unique kmer that diverges from this iteration's kmer
                group_tip_kmer = self.species[ next(group.__iter__()) ].tip_kmer()

                ### Assurance check:
                for i in group:
                    assert group_tip_kmer == self.species[i].tip_kmer()
                ### End check

                mut_kmer = self.kmer_generator.random_mutation(group_tip_kmer)
                for sid in group:
                    self.species[sid].extend(mut_kmer, diverge=True)
            else:
                for sid in group:
                    self.species[sid].extend(kmer)

        self.iteration += 1

    def save_fasta(self):
        """Saves the species sequences to a single fasta file"""
        pass

    def print_seqs(self):
        for s in self.species.values:
            print("> Species #{}\n  {}\n".format(s.sid, s.seq))

def main():
    km = KmerManager(9)
    sm = SpeciesManager(8, km)
    for i in range (45):
        if i % 5 == 0:
            print(sm.groups)
        sm.iterate()
    print("\n".join([i.seq for i in sm.species.values()]))
    #gg = km.generator.generate()
    #for i in range(10000):
    #    print(i, next(gg))

class TestSpeciesContainer(object):
    def test_modes(self):
        k = 3
        c = SpeciesContainer(1, k)
        c.extend("ATCG")
        c.extend("A")
        c.extend("G", diverge=True)
        assert c.cooldown_timer == k and c.mode == Mode.COOLDOWN
        c.extend("A")
        c.extend("A")
        c.extend("G", diverge=True)
        c.extend("A")
        c.extend("A")
        c.extend("A")
        assert c.cooldown_timer == 0 and c.mode == Mode.READY
        assert c.seq == "ATCGAGAAGAAA"

if __name__ == "__main__":
    km = KmerManager(9)
    sm = SpeciesManager(4, km)
    #main()
