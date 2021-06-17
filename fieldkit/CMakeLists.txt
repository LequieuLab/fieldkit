cmake_minimum_required(VERSION 3.10)

# copy/symlink python modules to the build directory to make it a working python package

file(GLOB files "*py") # option to glob

#set(files __init__.py
#          openfts.py
#          new_openfts.py
#          python_to_json.py
#   )
message(${files})

foreach(file ${files})
    # copy or symlink depending on cmake version
    if(${CMAKE_VERSION} VERSION_LESS "3.14.0") 
      file(COPY ${file}
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/)
    else()
      # symlink -- make life easier if user is real-time updating python files
      file(CREATE_LINK ${file}
           ${CMAKE_CURRENT_BINARY_DIR}/${file} SYMBOLIC)
    endif()
endforeach()

