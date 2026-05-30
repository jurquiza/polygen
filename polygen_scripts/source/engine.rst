engine module
=============

Tas and tigRNA additions
------------------------

PolyGEN supports TIGR-Tas guide RNA design through the TasA and TasH
architectures. TasA uses two programmable 9 nt spacers and TasH uses two
programmable 8 nt spacers. The mature tigRNA is built as edge repeat 5 prime,
spacer A, loop repeat, spacer B and edge repeat 3 prime.

The Tas page uses :func:`engine.tas_guide_design` with exact target windows for
multiplexing. Exact windows are 18 bp for TasA and 16 bp for TasH.
Exact-window entries may be separated by whitespace, commas, semicolons or
vertical bars.

Designed tigRNAs are passed into the assembly pipeline with the ``tigRNA``
polycistron type and input records of the form::

    tigRNA;sequence0|tigRNA;sequence1|tigRNA;sequence2

In ``tigRNA`` mode, :func:`engine.PTGbldr` treats every input sequence as a
mature tigRNA and does not add tRNAs or a Cpf1 direct repeat. During Golden Gate
optimization, internal overhangs are chosen from spacer B only for recognized
TasA/TasH scaffolds. Border oligos are target-specific and are always reported.
The melting temperature reported for tigRNA oligo pairs is calculated from the
shared forward/reverse overlap used for polymerase fill-in, not from the full
oligo sequence.

The GenBank output generated downstream of :func:`engine.scarless_gg` annotates
recognized tigRNAs with edge repeat 5 prime, spacer A, loop repeat, spacer B
and edge repeat 3 prime features.

.. automodule:: engine
    :members:
    :undoc-members:
    :show-inheritance:
