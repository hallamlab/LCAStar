Usage
=======

.. code-block:: python

    from lcastar import LcaStar, Lineage

Taxonomic inputs can be obtained from any source

Via scientific name (genus species)
-----------------------------------

Conversion taxonomic hierarchy is performed with the `ETE toolkit <https://github.com/etetoolkit/ete>`_.

.. code-block:: python

    orf_hits = [
        "Muribaculaceae bacterium",
        "Muribaculaceae bacterium",
        "Bacteroidales bacterium",
        "Muribaculaceae bacterium",
        "Alistipes senegalensis",
    ]

    tree = LcaStar()
    for sci_name in orf_hits:
        lin = Lineage.FromSciName(sci_name)
        assert lin is not None
        tree.NewObservation(lin)

Via NCBI taxonomy ID
--------------------

.. code-block:: python

    orf_hits = [
        2498093,
        2498093,
        2030927,
        2498093,
        1288121,
    ]

    tree = LcaStar()
    for tax_id in orf_hits:
        lin = Lineage.FromTaxID(tax_id)
        assert lin is not None
        tree.NewObservation(lin)

Via custom taxonomy
-------------------

Custom hierarchies can be used by providing traversal paths from the root node to each hit as a list of
tuples in the format `(rank, name)`. Paths need not be complete or reach the leaf. LCA\* will do its best
to construct the needed portions of the reference hierarchy. In this mode, LCA* is not limited to
biological taxonomies.

.. code-block:: python

    orf_hits = [
        [
            ("D", "Bacteria"),
            ("C", "Gammaproteobacteria"),
            ("G", "Escherichia"),
        ],
        [
            ("D", "Bacteria"),
            ("C", "Gammaproteobacteria"),
            ("G", "Pseudomonas"),
        ],
            [
            ("D", "Bacteria"),
            ("C", "Gammaproteobacteria"),
        ],
    ]

    tree = LcaStar()
    for path_from_root in orf_hits:
        lin = Lineage(path_from_root)
        tree.NewObservation(lin)

Example output
---------------

.. code-block:: python

    for node in tree.BestLineage():
        print(node.level, node.name, node.fraction_votes, node.p_value)

Expected output of the examples using NCBI taxonomy ID and scientific name.

.. code-block:: text

    superkingdom Bacteria 1.0 0.08273697918531309
    clade FCB group 1.0 0.08273697918531309
    clade Bacteroidota/Chlorobiota group 1.0 0.08273697918531309
    phylum Bacteroidota 1.0 0.08273697918531309
    class Bacteroidia 1.0 0.08273697918531309
    order Bacteroidales 1.0 0.08273697918531309
    species Bacteroidales bacterium 0.2 1.0
