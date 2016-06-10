#!/bin/bash

function makeACADOProject()

{
    projectBase="/home/arog/Desktop/"
    projectName=$1/
    projectDirectory=$projectBase$projectName

    if [ ! -d "$projectDirectory" ]
    then
	mkdir $projectDirectory
	echo "Directory $projectDirectory created"
    else
	echo "Directory $projectDirectory already exists"
    fi

    acadoRoot="/home/arog/Desktop/ACADOtoolkit"
    findAcadoFile="$acadoRoot/cmake/FindACADO.cmake"
    cp $findAcadoFile $projectDirectory

    baseCmakeFile="/home/arog/Desktop/ACADO_Aux/CMakeLists_BASE.txt"
    mf="CMakeLists.txt"
    cmakeFile=$projectDirectory$mf
    if [ -e "$cmakeFile" ]
    then
	rm $cmakeFile
	echo "$cmakeFile deleted"
    fi

    cp $baseCmakeFile $cmakeFile
    
	## Experimental
	pyFile="modifyCMakeFile.py"
	cp $pyFile $projectDirectory
	newPyFile=$projectDirectory$pyFile
	pyArg1="$2"
	pyArg2="$3"
	cd $projectDirectory
	python $pyFile $2 $3
	
	echo "Successfully built new project $1 at $projectDirectory, write CPP source files then compile using the script mkbuild.sh"
}

makeACADOProject rocketProblem rocket rocket.cpp
