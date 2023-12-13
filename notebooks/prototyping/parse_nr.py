from pathlib import Path
import re

# this is prolly not worth it...

def regex(r, s):
    for m in re.finditer(r, s):
        yield s[m.start():m.end()]

NR = Path("/home/tony/workspace/resources/nr/nr.faa")
OUT = Path("./cache/nr_meta.tsv")

with open(OUT, "w") as out:
    out.write("id\tname\ttax\n")
    with open(NR) as f:
        i = 0
        for l in f:
            if i % 100_000 == 0: print(i, end="\r")
            if not l.startswith(">"): continue
            i += 1
            
            try:
                if "[" in l and "]" in l:
                    id, rest = l[1:-1].split(" ", maxsplit=1)
                    name, rest = rest.split("[", maxsplit=1)
                    tax = " ".join(rest.split(" ")[:2])
                    tax = tax.replace("[", "").replace("]", "")
                    # tax = next(regex(r"([^\s]+\s?){1,2}", rest))
                    # tax = ''.join(c for c in tax if c not in "[]")
                else:
                    continue
                
                out.write(f"{id}\t{name}\t{tax}\n")
            except (StopIteration, ValueError) as e:
                print(i, e, l)
                break
