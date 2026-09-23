# Simple build for the JAGS Jones-Pewsey module.
#
# Usage:
#   make
#   sudo make install
#
# Override JAGS paths if needed:
#   make JAGS_INCLUDEDIR=/usr/local/include JAGS_LIBDIR=/usr/local/lib

JAGS_INCLUDEDIR ?= /usr/include/JAGS
JAGS_LIBDIR ?= /usr/lib/x86_64-linux-gnu
JAGS_MAJOR ?= 4

CXX ?= g++
CXXFLAGS = -std=c++11 -fPIC -O2 -Wall -I$(JAGS_INCLUDEDIR) -Isrc
LDFLAGS = -shared
LDLIBS = -L$(JAGS_LIBDIR) -ljags -ljrmath

SOURCES = src/jonespewsey.cc \
          src/distributions/DJonesPewsey.cc \
          src/functions/DJonesPewseyFun.cc
TARGET = jonespewsey.so
MODDIR = $(JAGS_LIBDIR)/JAGS/modules-$(JAGS_MAJOR)

.PHONY: all clean install test

all: $(TARGET)

$(TARGET): $(SOURCES) src/JonesPewseyMath.h src/distributions/DJonesPewsey.h \
           src/functions/DJonesPewseyFun.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SOURCES) $(LDLIBS)

install: $(TARGET)
	install -d $(DESTDIR)$(MODDIR)
	install -m 644 $(TARGET) $(DESTDIR)$(MODDIR)/

test: $(TARGET)
	python3 tests/test_jonespewsey.py

clean:
	rm -f $(TARGET)
