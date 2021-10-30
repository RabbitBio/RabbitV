DIR_INC := ./inc
DIR_SRC := ./src
DIR_OBJ := ./obj
DIR_IO := ./src/io

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin
INCLUDE_DIRS ?=
LIBRARY_DIRS ?=

SRC := $(wildcard ${DIR_SRC}/*.cpp)
IO := $(wildcard ${DIR_IO}/*.cpp)
OBJ := $(patsubst %.cpp,${DIR_OBJ}/%.o,$(notdir ${SRC}))
OBJ += $(patsubst %.cpp,${DIR_IO}/%.o,$(notdir ${IO}))

TARGET := RabbitV

BIN_TARGET := ${TARGET}

CXX ?= g++
CXXFLAGS := -std=c++11 -g -O3 -flto -funroll-loops -mavx512vl -mavx512bitalg -mavx512f -mavx512bw -mavx512vbmi2 -I${DIR_INC} $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) ${CXXFLAGS}
LIBS := -lz -lpthread
LD_FLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS) $(LD_FLAGS)


${BIN_TARGET}:${OBJ}
	$(CXX) -flto -g $(OBJ) -o $@ $(LD_FLAGS)

#${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp make_obj_dir
#	$(CXX) -c $< -o $@ $(CXXFLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.cpp 
	$(CXX) -c $< -o $@ $(CXXFLAGS)

.PHONY:clean
clean:
	@if test -d $(DIR_OBJ) ; \
	then \
		find $(DIR_OBJ) -name *.o -delete; \
		find $(DIR_IO) -name *.o -delete; \
	fi
	@if test -e $(TARGET) ; \
	then \
		rm $(TARGET) ; \
	fi

make_obj_dir:
	@if test ! -d $(DIR_OBJ) ; \
	then \
		mkdir $(DIR_OBJ) ; \
	fi

install:
	install $(TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."
