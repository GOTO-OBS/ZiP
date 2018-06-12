#import catsHTM2
from wsdb import WSDB
from astropy.table import Table

def ctsrch(RA, DEC):
 client = WSDB(database='wsdb', host="goto-observatory.org", user="goto", password="A PASSWORD")

#print(client.catalogues)
 return(client.radial_query(catalogue_name="apassdr9_main", ra=RA, dec=DEC, radius=0.05))

#RA= 45.212
#DEC= -34.113
#LIST = ctsrch(RA,DEC)

#print(LIST)

#for i in range(0, len(LIST)):
# print(type(float(LIST[i][2])))
# print(LIST[i][4])
# print(LIST[i][10]) #gmag
