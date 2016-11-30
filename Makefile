#Copyright 2016 Dmitry N. Razdoburdin.

#This file is part of TG_3D. TG_3D is a program that calculats transient growth of linear perturbations in 3D shearing box by power iterations.
#TG_3D is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. TG_3D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with TG_3D; if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

SOURSE_DIR=./sourses
CONF_DIR=./configs
OBJECTS_DIR=./objects

include $(CONF_DIR)/params.conf
include $(CONF_DIR)/compile_keys.conf
include $(CONF_DIR)/link_keys.conf

ifeq ($(BACKGROUND),isothermal)
	BOUNDARY=infinite
endif

DEFINES=
ifeq ($(SIGNAL2),yes)
	DEFINES+=-D SIGNAL2
endif
ifeq ($(SIGNAL3),yes)
	DEFINES+=-D SIGNAL3
endif
ifeq ($(SIGNAL4),yes)
	DEFINES+=-D SIGNAL4
endif
ifeq ($(TEST_OF_CONJUGATION),yes)
	DEFINES +=-D TEST_OF_CONJUGATION
endif
ifeq ($(G_OUTPUT),full)
	DEFINES +=-D G_OUTPUT_FULL
endif
ifeq ($(LOG_OUTPUT),stderr)
	DEFINES +=-D LOG=stderr
else
	DEFINES +=-D LOG=data.LogFile
	DEFINES +=-D LOGFILENAME=$(LOG_OUTPUT)
endif

DEFINES_OMP=
ifeq ($(VECTORIZE),yes)
	DEFINES_OMP+=-D SIMD=simd
endif
ifeq ($(VECTORIZE),no)
	DEFINES_OMP+=-D SIMD=
endif
ifeq ($(BOUNDARY),first)
	DEFINES_OMP+= -D WAIT=nowait
endif
ifeq ($(BACKGROUND),isothermal)
	DEFINES_OMP+= -D WAIT=nowait
endif
ifeq ($(BOUNDARY),second)
	DEFINES_OMP+= -D WAIT=
endif
ifeq ($(BOUNDARY),periodic)
	DEFINES_OMP+= -D WAIT=
endif

TG_3D : $(OBJECTS_DIR)/main.o $(OBJECTS_DIR)/methods.o $(OBJECTS_DIR)/procedures.o $(OBJECTS_DIR)/functions.o $(OBJECTS_DIR)/functions_background.o  $(OBJECTS_DIR)/methods_metric.o $(OBJECTS_DIR)/methods_boundary.o $(CONF_DIR)/compile_keys.conf $(CONF_DIR)/link_keys.conf Makefile
	$(CXX) $(LDFLAGS) $(OBJECTS_DIR)/main.o $(OBJECTS_DIR)/methods.o $(OBJECTS_DIR)/procedures.o $(OBJECTS_DIR)/functions.o $(OBJECTS_DIR)/functions_background.o $(OBJECTS_DIR)/methods_metric.o $(OBJECTS_DIR)/methods_boundary.o -o TG_3D $(LDLIBS)

$(OBJECTS_DIR)/main.o : $(SOURSE_DIR)/main.cpp $(SOURSE_DIR)/classes.h $(SOURSE_DIR)/procedures.h $(CONF_DIR)/compile_keys.conf
	mkdir -p $(OBJECTS_DIR)
	$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/main.cpp -o $(OBJECTS_DIR)/main.o

$(OBJECTS_DIR)/methods.o : $(SOURSE_DIR)/methods.cpp $(SOURSE_DIR)/classes.h $(SOURSE_DIR)/functions.h $(CONF_DIR)/compile_keys.conf $(CONF_DIR)/link_keys.conf
	$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/methods.cpp -o $(OBJECTS_DIR)/methods.o $(LDLIBS) $(DEFINES) $(DEFINES_OMP)

$(OBJECTS_DIR)/functions.o : $(SOURSE_DIR)/functions.cpp $(SOURSE_DIR)/classes.h $(CONF_DIR)/compile_keys.conf
	$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/functions.cpp -o $(OBJECTS_DIR)/functions.o

$(OBJECTS_DIR)/functions_background.o : $(SOURSE_DIR)/functions_$(BACKGROUND).cpp $(CONF_DIR)/compile_keys.conf $(CONF_DIR)/link_keys.conf
	$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/functions_$(BACKGROUND).cpp -o $(OBJECTS_DIR)/functions_background.o

$(OBJECTS_DIR)/methods_metric.o : $(SOURSE_DIR)/methods_$(METRIC).cpp $(SOURSE_DIR)/classes.h $(SOURSE_DIR)/functions.h $(CONF_DIR)/compile_keys.conf $(CONF_DIR)/link_keys.conf
	$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/methods_$(METRIC).cpp -o $(OBJECTS_DIR)/methods_metric.o $(LDLIBS) $(DEFINES_OMP)

$(OBJECTS_DIR)/methods_boundary.o : $(SOURSE_DIR)/methods_$(BOUNDARY).cpp $(SOURSE_DIR)/classes.h $(SOURSE_DIR)/functions.h $(CONF_DIR)/compile_keys.conf $(CONF_DIR)/link_keys.conf
		$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/methods_$(BOUNDARY).cpp -o $(OBJECTS_DIR)/methods_boundary.o

$(OBJECTS_DIR)/procedures.o : $(SOURSE_DIR)/procedures.cpp $(SOURSE_DIR)/classes.h $(CONF_DIR)/compile_keys.conf $(CONF_DIR)/link_keys.conf
		$(CXX) -c $(CXXFLAGS) $(SOURSE_DIR)/procedures.cpp -o $(OBJECTS_DIR)/procedures.o $(LDLIBS)

recompile :
	rm -rf $(OBJECTS_DIR)/*.o TG_3D
	make TG_3D

task : TG_3D $(CONF_DIR)/params.conf
	./TG_3D $(KEYS)

profiling : TG_3D $(CONF_DIR)/params.conf
	rm -f callgrind.out.*
#	valgrind --leak-check=yes --track-origins=yes ./TG_3D $(KEYS)
#	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes -v ./TG_3D $(KEYS)
	valgrind --tool=callgrind ./TG_3D $(KEYS)
