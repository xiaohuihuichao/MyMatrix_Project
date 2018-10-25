CMAKE_MINIMUM_REQUIRED( VERSION 3.0 )



PROJECT( Project )




INCLUDE_DIRECTORIES( MyMatrix/include )


AUX_SOURCE_DIRECTORY( MyMatrix/src DIR_SRC )

AUX_SOURCE_DIRECTORY( MyMatrix DIR_MAIN )



ADD_SUBDIRECTORY( module/0.Linear_Regression )


FILE(GLOB_RECURSE HEADERS *.h *.hpp)
SOURCE_GROUP("Include" FILES ${HEADERS} ) 




ADD_EXECUTABLE( 00.MyMatrix  ${DIR_SRC} 
			${DIR_MAIN}
			${HEADERS} )
