# enforces C89, to ensure it will even compile on a Microsoft C compiler
CC             = gcc
CFLAGS         = -Wall -std=c89 -pedantic -O2 -Wshadow -Wpointer-arith \
	-Wcast-qual -Wcast-align -Wwrite-strings \
	-Wstrict-prototypes -Wmissing-prototypes -Wnested-externs -Wno-long-long \
	-DMAXSEQ=$(MAX_SEQ) -g #-pg
LDFLAGS        = -lm #-pg

# Customize Weeder program here
MAX_SEQ        = 100000
ADVISER        = adviser
WEEDER_TFBS    = weederTFBS
WEEDERLAUNCHER = weederlauncher
LOCATOR        = locator
#LAUNCHERFLAGS  = -DADVISER_CMD="\"./$(ADVISER)\"" -DWEEDER_TFBS_CMD="\"./$(WEEDER_TFBS)\""
LAUNCHERFLAGS  = -DADVISER_CMD="\"$(ADVISER)\"" -DWEEDER_TFBS_CMD="\"$(WEEDER_TFBS)\""
INSTALL_DIR = /Users/wwu/local/bin

all: $(WEEDER_TFBS) $(ADVISER) $(LOCATOR) $(WEEDERLAUNCHER)

install: all
	cp $(WEEDER_TFBS) $(ADVISER) $(LOCATOR) $(WEEDERLAUNCHER) $(INSTALL_DIR)

clean:
	rm -f $(WEEDER_TFBS) $(ADVISER) $(LOCATOR) $(WEEDERLAUNCHER) \
		*.o weeder_util_test

test: weeder_util_test
	./weeder_util_test

weeder_util_test: weeder_util_test.o CuTest.o weeder_util.o
	$(CC) weeder_util_test.o weeder_util.o CuTest.o -o weeder_util_test $(LDFLAGS)

weeder_util_test.o: weeder_util_test.c weeder_util.h
	$(CC) weeder_util_test.c -c

$(WEEDER_TFBS): weeder_util.o weederTFBS.o
	$(CC) weederTFBS.o weeder_util.o -o $(WEEDER_TFBS) $(LDFLAGS)

$(ADVISER): weeder_util.o adviser.o
	$(CC) adviser.o weeder_util.o -o $(ADVISER) $(LDFLAGS)

$(LOCATOR): weeder_util.o locator.o
	$(CC) locator.o weeder_util.o -o $(LOCATOR) $(LDFLAGS)

$(WEEDERLAUNCHER): weederlauncher.o weeder_util.o
	$(CC) weederlauncher.o weeder_util.o -o $(WEEDERLAUNCHER) $(LDFLAGS)

weederlauncher.o: weeder_util.h weederlauncher.c
	$(CC) weederlauncher.c -c $(CFLAGS) $(LAUNCHERFLAGS)

weeder_util.o: weeder_util.h weeder_util.c
	$(CC) weeder_util.c -c $(CFLAGS)

weederTFBS.o: weeder_util.h weederTFBS.c
	$(CC) weederTFBS.c -c $(CFLAGS)

locator.o: weeder_util.h locator.c
	$(CC) locator.c -c $(CFLAGS)

adviser.o: weeder_util.h adviser.c
	$(CC) adviser.c -c $(CFLAGS)
