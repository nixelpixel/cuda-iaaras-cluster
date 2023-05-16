CC = mpicc
SRCPATH = ./src
OUTPATH = ./build
MPILIBS = -lmpi -lm
MATHLIBFLAGS = -D _USE_MATH_DEFINES
CFLAGS = -Wall -Wextra -Wpedantic -std=c99 -g

CC_OBJ = $(CC) -c $(CFLAGS) $(SRCPATH)
TO = -o $(OUTPATH)


all: lib_cluster.o lib_blocklist.o ctldev.o corrs.o main.o
	$(CC) $(wildcard $(OUTPATH)/*.o) $(MPILIBS) $(FFTW3LIBS) -o $(OUTPATH)/cluster

main.o: $(SRCPATH)/cluster_main.c $(SRCPATH)/cluster_lib.h
	$(CC_OBJ)/cluster_main.c $(MPILIBS) $(TO)/main.o

lib_cluster.o: $(SRCPATH)/cluster_lib.c $(SRCPATH)/cluster_lib.h
	$(CC_OBJ)/cluster_lib.c $(MPILIBS) $(TO)/lib_cluster.o

lib_blocklist.o: $(SRCPATH)/blocklist_lib.c $(SRCPATH)/blocklist_lib.h $(SRCPATH)/cluster_lib.h
	$(CC_OBJ)/blocklist_lib.c $(MPILIBS) $(TO)/lib_blocklist.o

ctldev.o: $(SRCPATH)/ctldev_main.c $(SRCPATH)/blocklist_lib.h
	$(CC_OBJ)/ctldev_main.c $(MPILIBS) -l $(OUTPATH)/lib_CTF.o $(TO)/ctldev.o

corrs.o: lib_CTF $(SRCPATH)/corr_main.c $(SRCPATH)/cluster_lib.h
	$(CC_OBJ)/corr_main.c $(MPILIBS) -l $(OUTPATH)/lib_CTF.o $(TO)/corrs.o

lib_CTF: $(SRCPATH)/CFT_lib.c $(SRCPATH)/CFT_lib.h
	$(CC_OBJ)/CFT_lib.c $(MATHLIBFLAGS) $(TO)/lib_CTF.o
	ar rcs $(OUTPATH)/lib_CTF.a $(OUTPATH)/lib_CTF.o


clean:
	rm $(wildcard $(OUTPATH)/*.o) $(wildcard $(OUTPATH)/*.exe)


#all: lib_cluster lib_blocklist ctldev corr
#	$(CC) -c $(CFLAGS) $(wildcard $(OUTPATH)/*.obj) -o $(OUTPATH)/out.exe
#lib_cluster:
#	$(CC) -c $(CFLAGS) $(SRCPATH)/cluster_lib.c -o $(OUTPATH)/cluster_lib.obj
#lib_blocklist:
#	$(CC) -c $(CFLAGS) $(SRCPATH)/blocklist_lib.c -o $(OUTPATH)/blocklist_lib.obj
#ctldev:
#	$(CC) -c $(CFLAGS) $(SRCPATH)/ctldev_main.c -o $(OUTPATH)/ctldev.obj
#corr:
#	$(CC) -c $(CFLAGS) $(SRCPATH)/corr_main.c -o $(OUTPATH)/corr.obj
