# This script identifies text file in a given directory and appends a new line at the end of each file. Developed to keep MDOODZ .txt files clean such that all parameters are read in properly.

function cleanTxtFiles()
    # Set vars
    targetDir = "/Users/lcandiot/Developer/Julia/BasicScripts/data/"
    suffix    = "txt"
    strEOF    = "/********** END OF FILE **********/"
    # Get directory list
    for f in filter(x->endswith(x,suffix),readdir(targetDir))
        # Open file just to read the existing lines
        file = open("$(targetDir)$(f)") 
        lines = readlines(file)
        close(file)
        # Append EOF string in case it is missing
        if lines[end]!=strEOF
            push!(lines, strEOF)
        end
        # Write new version of the .txt file
        file = open("$(targetDir)$(f)_copy", "w") 
        for line in lines
            println(file, line)
        end
        close(file)
        # Remove old and update file names
        run(`rm -r $(targetDir)$(f)`)
        run(`mv $(targetDir)$(f)_copy $(targetDir)$(f)`)
    end
end

cleanTxtFiles()