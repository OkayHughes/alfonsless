function try_exit(st)
    try
        run(st)
    catch
	"try_exit encountered error"
        exit(1)
    end
    exit
end
