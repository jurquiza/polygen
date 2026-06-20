engine module
=============

Tas and tigRNA additions
------------------------

PolyGEN supports TIGR-Tas guide RNA design through the TasA and TasH
architectures. TasA uses two programmable 9 nt spacers and TasH multiplexing
uses two programmable 8 nt spacers. The mature tigRNA is built as edge repeat 5 prime,
spacer A, loop repeat, spacer B and edge repeat 3 prime.
The default TasA scaffold is AACCG, spacer A, AGTAACCCC, spacer B, AGTG.
The default TasH scaffold is GAGCGAGTTA, spacer A, AAAACAATCA,
spacer B, AAGCGAGCCA.

The Tas page uses :func:`engine.tas_guide_design` with exact target windows for
multiplexing. Exact windows are 18 bp for TasA and 16 bp for TasH. For TasH
multiplexing, the first 8 bp are embedded as spacer A and the reverse
complement of the second 8 bp is embedded as spacer B.
Exact-window entries may be separated by whitespace, commas, semicolons or
vertical bars.

For adjacent asymmetric Tas units, :func:`engine.PTGbldr` tracks the upstream
right edge plus the downstream left edge as required edge-junction context. For
the current defaults this junction is AGTG-AACCG in TasA arrays and
AAGCGAGCCA-GAGCGAGTTA in TasH arrays.

Designed tigRNAs are passed into the assembly pipeline with the ``tigRNA``
polycistron type and input records of the form::

    tigRNA;sequence0|tigRNA;sequence1|tigRNA;sequence2

In ``tigRNA`` mode, :func:`engine.PTGbldr` treats every input sequence as a
mature tigRNA and does not add tRNAs or a Cpf1 direct repeat. During Golden Gate
optimization, internal overhangs are chosen from spacer A or spacer B for
recognized TasA/TasH scaffolds while avoiding fixed edge and loop sequences.
Border oligos are target-specific and are always reported.
Oligo arms are dynamically balanced up to 60 bp when needed so short TasA/TasH
fragments keep enough shared overlap for assembly. The melting temperature
reported for tigRNA oligo pairs is calculated from the shared forward/reverse
overlap used for polymerase fill-in, not from the full oligo sequence.

The GenBank output generated downstream of :func:`engine.scarless_gg` annotates
recognized tigRNAs with edge repeat 5 prime, spacer A, loop repeat, spacer B
and edge repeat 3 prime features. Primer-binding annotations include the 4 bp
Golden Gate overhang when it matches the final construct, but exclude the Type
IIS recognition site.

.. automodule:: engine
    :members:
    :undoc-members:
    :show-inheritance:
