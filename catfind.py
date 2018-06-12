#import catsHTM2
from wsdb import WSDB
from astropy.table import Table

def ctsrch(RA, DEC):
 client = WSDB(database='wsdb', host="goto-observatory.org", user="goto", password="A PASSWORD")

#print(client.catalogues)
 return(client.radial_query(catalogue_name="apassdr9_main", ra=RA, dec=DEC, radius=0.05))
