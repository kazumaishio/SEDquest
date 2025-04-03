import SEDquest

obj = SEDquest.SEDgenerator('a')

from astropy.table import QTable
testtable = QTable.read('./collarea_magic.ecsv')
testtable