LIBDIR = $(DESTDIR)/usr/share/pyshared
BINDIR = $(DESTDIR)/usr/bin
CONFDIR = $(DESTDIR)/var/cmonkey2

clean:
	rm -f *.pyc[co] */*.py[co]

install:
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)
	mkdir -p $(CONFDIR)
	cp -r cmonkey $(LIBDIR)/
	cp cmonkey.py $(BINDIR)/cmonkey2
	cp config/* $(CONFDIR)/

uninstall:
	rm -rf $(LIBDIR)/cmonkey2
	rm -f $(BINDIR)/cmonkey2
	rm -rf $(CONFDIR)

