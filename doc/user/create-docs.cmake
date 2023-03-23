# This script generates the documentation for the Crux commands using the
# 'crux create-docs' command.

# Create the list of documents to create
execute_process(
  COMMAND ${CRUX_PATH} create-docs --no-analytics T list
  RESULT_VARIABLE status
  OUTPUT_VARIABLE doc_list 
  ERROR_VARIABLE error_message
)
if (NOT ${status} EQUAL 0)
  message(
    FATAL_ERROR
    "Unable to create list of documents " 
    "\nError message: ${error_message}"
    "\nCreation of documents failed."
  )
endif (NOT ${status} EQUAL 0)

# Documents are separated by newlines, convert to CMake list
STRING(REGEX REPLACE "\n" ";" doc_list "${doc_list}")

# Generate each document
foreach (doc ${doc_list})
  execute_process(
    COMMAND ${CRUX_PATH} create-docs --no-analytics T ${doc}
      RESULT_VARIABLE status
      ERROR_VARIABLE error_message
      OUTPUT_FILE "${doc}.html"
  )
  if (NOT ${status} EQUAL 0)
    message(
      FATAL_ERROR
      "Unable to generate ${doc}.html" 
      "\nError message: ${error_message}"
      "\nCreation of documents failed."
    )
  endif (NOT ${status} EQUAL 0)
  message(STATUS "Created ${doc}.html")
  
  # Copy file to source doc/commands directory.
  message(STATUS "Copying doc files back to source doc/command directory: ${DOC_DIR}")
  execute_process(
    COMMAND ${CMAKE_COMMAND} 
      -E copy ${doc}.html ${DOC_DIR}
      RESULT_VARIABLE status
      ERROR_VARIABLE error_message
  )
  if (NOT ${status} EQUAL 0)
    message(
      FATAL_ERROR
      "Unable to copy ${doc}.html to ${DOC_DIR}" 
      "\nError message: ${error_message}"
      "\nCreation of documents failed."
    )
  endif (NOT ${status} EQUAL 0)
endforeach (doc ${doc_list})

# Create the parameter table
execute_process(
  COMMAND ${CRUX_PATH} create-docs --no-analytics T param-table
    RESULT_VARIABLE status
    ERROR_VARIABLE error_message
    OUTPUT_FILE param-table.html
)
if (NOT ${status} EQUAL 0)
  message(
    FATAL_ERROR
    "Unable to generate param-table.html" 
    "\nError message: ${error_message}"
    "\nCreation of param-table failed."
  )
endif (NOT ${status} EQUAL 0)
message(STATUS "Created param-table.html")

# Copy file to source doc/commands directory.
execute_process(
  COMMAND ${CMAKE_COMMAND} 
    -E copy param-table.html ${DOC_DIR}
    RESULT_VARIABLE status
    ERROR_VARIABLE error_message
)
if (NOT ${status} EQUAL 0)
  message(
    FATAL_ERROR
    "Unable to copy param-table.html to ${DOC_DIR}" 
    "\nError message: ${error_message}"
    "\nCreation of param-table failed."
  )
endif (NOT ${status} EQUAL 0)

# Create the default.params file
execute_process(
  COMMAND ${CRUX_PATH} create-docs --no-analytics T default-params
    RESULT_VARIABLE status
    ERROR_VARIABLE error_message
    OUTPUT_FILE default.params
)
if (NOT ${status} EQUAL 0)
  message(
    FATAL_ERROR
    "Unable to generate default.params"
    "\nError message: ${error_message}"
    "\nCreation of default.params failed."
  )
endif (NOT ${status} EQUAL 0)
message(STATUS "Created default.params")

# Copy file back to source doc/commands directory.
execute_process(
  COMMAND ${CMAKE_COMMAND} 
    -E copy default.params ${DOC_DIR}
    RESULT_VARIABLE status
    ERROR_VARIABLE error_message
)
if (NOT ${status} EQUAL 0)
  message(
    FATAL_ERROR
    "Unable to copy default.params to ${DOC_DIR}" 
    "\nError message: ${error_message}"
    "\nCreation of param-table failed."
  )
endif (NOT ${status} EQUAL 0)
