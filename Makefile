LIBDIR = $(DESTDIR)/usr/share/pyshared
BINDIR = $(DESTDIR)/usr/bin
CONFDIR = $(DESTDIR)/var/cmonkey-python

clean:
	rm -f *.pyc[co] */*.py[co]

install:
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)
	mkdir -p $(CONFDIR)
	cp -r cmonkey $(LIBDIR)/
	cp cmonkey.py $(BINDIR)/cmonkey-python
	cp config/* $(CONFDIR)/

uninstall:
	rm -rf $(LIBDIR)/cmonkey
	rm -f $(BINDIR)/cmonkey-python
	rm -rf $(CONFDIR)
