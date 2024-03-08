from __future__ import annotations
from dataclasses import dataclass, field
from typing import Iterable
import numpy as np
from scipy.stats import chi2
# http://etetoolkit.org/download/
# http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html#setting-up-a-local-copy-of-the-ncbi-taxonomy-database
# warning: ete tk inserts itself into ~/.local/share/ete/
from ete4 import NCBITaxa

_ete_ncbi = None
def InitNcbi() -> NCBITaxa:
    global _ete_ncbi
    if _ete_ncbi is None:
        _ete_ncbi = NCBITaxa()
    return _ete_ncbi

class Lineage:
    @classmethod
    def FromTaxID(cls, taxid: int) -> Lineage|None:
        ncbi = InitNcbi()
        id_lin = ncbi.get_lineage(taxid)
        if id_lin is None: return None
        names = ncbi.get_taxid_translator(id_lin)
        ranks = ncbi.get_rank(id_lin)
        return cls((names[id], ranks[id]) for id in id_lin if ranks[id] != "no rank")

    @classmethod
    def FromSciName(cls, sci_name: str) -> Lineage|None:
        ncbi = InitNcbi()
        _get = lambda x: next(iter(x.values()))
        tax_data = ncbi.get_name_translator([sci_name])
        if len(tax_data)==0: # retry with genus
            tax_data = ncbi.get_name_translator([sci_name.split(" ")[0]])
        if len(tax_data)==0:
            return None
        else:
            taxid = _get(tax_data)[0]
        return cls.FromTaxID(taxid)

    def __init__(self, lineage: Iterable[tuple[str, str]]|None = None) -> None:
        """@lineage is iterable of (name, level) tuples"""
        self._i = 0
        self._lineage: list[tuple[str, str]] = [] if lineage is None else list(lineage)

    def Add(self, name: str, level: str):
        self._lineage.append((name, level))

    def __iter__(self):
        self._i = 0
        return self
    
    def __next__(self):
        if self._i >= len(self._lineage): raise StopIteration
        self._i += 1
        return self._lineage[self._i-1]
    
@dataclass
class ResultNode:
    name: str
    level: str
    entropy: float
    cumulative_votes: int
    fraction_votes: float
    p_value: float

class LcaStar:
    @dataclass
    class Node:
        name: str
        level: str
        entropy: float =    field(default_factory=lambda: 0.0)
        count: int =        field(default_factory=lambda: 0)
        sum_entropy: float =field(default_factory=lambda: 0.0)
        sum_count: int =    field(default_factory=lambda: 0)
        parent: LcaStar.Node|None = field(default_factory=lambda: None)

        def __str__(self) -> str:
            return f"<{self.name}, entropy={self.sum_entropy}, count={self.sum_count}>"
        
        def __repr__(self) -> str:
            return self.__str__()
        
    def __init__(self) -> None:
        self.root = LcaStar.Node("root", "root")
        self.nodes: dict[str, LcaStar.Node] = {self.root.name:self.root}
        self.sum = 0
        self.counts = set()

    def AddPlacementsAtRoot(self, count: int):
        """decrease confidence to take into account, for example, unannotated ORFs on the contig"""
        self.root.count += count

    def NewObservation(self, lineage: Lineage):
        self.sum += 1
        parent = self.root
        for k, level in lineage:
            if k not in self.nodes:
                self.nodes[k] = LcaStar.Node(k, level, parent=parent)
            else:
                assert self.nodes[k].parent == parent, f"invalid lineage: the parent of [{k}] was previously said to be [{self.nodes[k].parent}], now [{parent.name}]"
            parent = self.nodes[k]
        leaf = parent # last node assigned in the loop
        leaf.count += 1

    def BestLineage(self) -> list[ResultNode]:
        # -----------------------------------------------------
        # update entropy

        for node in self.nodes.values():
            node.entropy = 0.0

        for node in self.nodes.values():
            if node.count == 0: continue

            self.counts.add(node.count)
            r_of_x = node.count/self.sum
            c = node.count
            entropy = r_of_x*np.log10(r_of_x)
            while node is not None:
                node.sum_entropy += entropy
                node.sum_count += c
                node = node.parent

        # -----------------------------------------------------
        # prepare for finding best lineage

        all_children: dict[str, list[LcaStar.Node]] = {}
        for node in self.nodes.values():
            if node.parent is None: continue
            if node.parent.name not in all_children:
                all_children[node.parent.name] = []
            all_children[node.parent.name].append(node)
        
        # for calculating p-value
        other_roots: list[LcaStar.Node] = []
        _last_len, _last_largest_other = 0, 0
        others: list[LcaStar.Node] = [] # going down tree, in case some counts are not at leaf
        def _find_second_largest_count():
            nonlocal _last_len, _last_largest_other
            _get_max = lambda: max(max(n.count for n in others), _last_largest_other)
            if _last_len == len(other_roots): return _get_max()
            todo = other_roots.copy()
            while len(todo)>0:
                node = todo.pop()
                if node.count > _last_largest_other:
                    _last_largest_other = node.count
                todo += all_children.get(node.name, [])
            _last_len = len(other_roots)
            return _get_max()

        # p-value is based on supremacy, where null hypothesis is:
        # there actially exists another taxon with a count larger than @node.count
        def _pvalue(node: LcaStar.Node, second_largest_count: int) -> float:
            M = float(second_largest_count) # "M" from supplementary
            X_k = float(node.sum_count) # "X_k" from supplementary
            if X_k <= M:
                t = 0.0 # trivial case
            else:
                first = 0
                if M > 0:
                    first = M * np.log10( (2 * M) / (M + X_k) )
                second = X_k * np.log10( (2 * X_k) / (M + X_k) )
                t = 2.0*(first + second)
            # "...it can also be shown that an approximate likelihood-ratio test at some significance level α
            # can be approximated by a χ-square distribution of degree one"
            # https://doi.org/10.1198/jasa.2009.tm08213
            degrees_of_freedom = 1
            return float(1 - chi2.cdf(t, degrees_of_freedom))

        # -----------------------------------------------------
        # find best lineage

        lineage: list[ResultNode] = []
        if self.root.name not in all_children: return lineage
        
        children = all_children[self.root.name]
        while True:
            entropy_ordered = sorted(children, key=lambda x: x.sum_entropy)
            best = entropy_ordered[0]
            if len(entropy_ordered)>1:
                for other in entropy_ordered[1:]:
                    other_roots.append(other)
            if best.parent is not None: others.append(best.parent) 

            lineage.append(ResultNode(
                name                = best.name,
                level               = best.level,
                entropy             = best.sum_entropy,
                cumulative_votes    = best.sum_count,
                fraction_votes      = best.sum_count/self.sum,
                p_value             = _pvalue(best, _find_second_largest_count()),
            ))
            if best.name not in all_children: break
            children = all_children[best.name]
        return lineage
