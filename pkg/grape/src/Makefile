CC="cc"
CFLAGS="-O"

all:
	$(CC) ${LDFLAGS} $(CFLAGS) -o ../bin/drcanon3  drcanon3.c
	$(CC) ${LDFLAGS} $(CFLAGS) -o ../bin/drtogap3  drtogap3.c
	$(CC) ${LDFLAGS} $(CFLAGS) -o ../bin/gap3todr  gap3todr.c

clean:
	rm -f ../bin/drcanon3 \
	      ../bin/drtogap3 \
              ../bin/gap3todr
