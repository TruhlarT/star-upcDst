cmake_minimum_required (VERSION 3.0)
project (upc-reader)

#sources for LIB_READER
set (SRCS
  src/Util.cxx
  src/EmbedMaker.cxx
  src/RecTree.cxx
  src/SaveRecTree.cxx
  src/Ana.cxx
  src/MainAna.cxx
  src/VertexStudy.cxx
  src/TofQA.cxx
  src/TofEff.cxx
  src/TrigEff.cxx 
  src/FullZB.cxx 
  src/ElasticAna.cxx
  src/RpMCAna.cxx
  src/Embedding.cxx
  src/PlotVertexStudy.cxx
  src/PlotUtil.cxx 
  src/PlotMainAna.cxx
  src/PlotMissingPt.cxx
  src/PlotBackgroundStudy.cxx 
  src/PlotGraniitti.cxx 
  src/PlotSysStudy.cxx 
  src/PlotGetEfficiency.cxx       
  src/PlotEmbedding.cxx
  src/PlotTrigEff.cxx 
  src/PlotTofQA.cxx  
  src/PlotProbOfRetainEvent.cxx
  src/PlotRPMCPlots.cxx 
  src/PlotPID.cxx 
  )





#name of LIB_READER
set (LIB_READER upc-reader)

#binary executables
set (BIN AnalysisManager PlotManager)

#UPC picoDsts library
set (UPC_LIB "${CMAKE_CURRENT_SOURCE_DIR}/../build/libstar-upc.so")
set (UPC_INC "${CMAKE_CURRENT_SOURCE_DIR}/../include")

#general compiler options
set (CMAKE_CXX_COMPILER /usr/bin/g++)
set (CMAKE_CXX_STANDARD 11)

#option -Wl,--no-as-needed is necessary only with root5
set (CMAKE_CXX_FLAGS "-Wall -g -Wl,--no-as-needed")

include_directories (include ${UPC_INC})



#ROOT section

#compiler flags
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_FLAGS_CMD)
string(REPLACE "\n" "" ROOT_FLAGS "${ROOT_FLAGS_CMD}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_FLAGS}")

#libraries
set(ROOT_DEPS Core EG Hist Physics RIO Tree MathCore RooFit RooFitCore Minuit) 

# Add the ROOT libraries
set(ROOT_LIBRARIES ${ROOT_LIBRARIES} -lRooFit -lRooFitCore -lMinuit)


#include dir
execute_process(COMMAND root-config --libdir OUTPUT_VARIABLE ROOT_LIBDIR_CMD)
string(REPLACE "\n" "" ROOT_LIBDIR "${ROOT_LIBDIR_CMD}")
link_directories (${ROOT_LIBDIR})

#end of ROOT section

#create library LIB_READER
add_library (${LIB_READER} SHARED ${SRCS})
target_link_libraries(${LIB_READER} ${ROOT_DEPS} )

#create executables in BIN
foreach(IBIN ${BIN})
  add_executable(${IBIN} src/${IBIN}.cxx)
  target_link_libraries(${IBIN} ${LIB_READER} ${UPC_LIB} ${ROOT_DEPS} )
endforeach(IBIN ${BIN})






















