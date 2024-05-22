# TODO

Random observations, things that could be fixed.

## my-dynamic-dropdown is slow

This is a substring search over a python list.  SQLite can do this
instead (LIKE '%...%'), and if the substring search is replaced witha
prefix search (LIKE '...%'), it should be a lot faster.

## download\_table\_TPMS\_without\_means is strange

It uses the global variable 'tissues' and ignores the tissue selection
in the UI.  This is probably unintentional. 

The code is also very convoluted and could probably be replaced with a
single SQL statement.


## zooming doesn't work

If a region is entered and "Update Transcripts" clicked, the diagram
isn't updated to the region.  Instead, transcripts overlapping the
region are queried, and the displayed region is the smallest one that
includes all of them.

If this was changed to display an arbitrary region, the diagramming code
would need to cope with transcripts that don't fit entirely into the
diagram.


## only one strand is shown

By default, only transcripts on one strand are shown, and the display is
flipped so transcripts run left-to-right.  It might make sense to show
antisense transcripts, too, but the diagram should clearly indicate
theyr direction.


## database is constantly reconnected

Should check if holding one global connection works.  It may be necessary to
set 'autocommit' or sprinkle explicit commits somewhere.  While at it, creating
the cursor objects is usually pointless.


## convoluted code that could be SQL

- calculate\_inner\_exons
- download\_table\_TPMS\_without\_means
- get\_TPM\_from\_tissues
- get\_TPM\_from\_tissues\_over\_transcripts


## get\_transcripts\_from\_gene needs refactoring

Depending on how it was triggered, this monolith of a function does six
similar, but subtly different things.  Due to copy-and-paste coding,
some of the differences are probably unintentional.  They are also
unexpected from a user perspective.

Commonalities should be factored out, redundancies removed, and
similarity be enforced through code reuse.

