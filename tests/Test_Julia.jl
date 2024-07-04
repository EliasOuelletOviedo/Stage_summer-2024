

filepath = "/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-18_TestSoftware/test save online-001/Bliq VMS/C1/test save online-001-C1-T001.tif"
filepath = "/Volumes/home/2024-07-02_real-time_VMS/test_VMS-001/Bliq VMS/C1/test_VMS-001-C1-T012.tif"

println(String(read(filepath, 3680)))

framerate = match(r"<M K=\"FrameRate\">(.*?)</M>", String(read(filepath, 3680))).captures[1]
sizex = match(r" SizeX=\"(.*?)\"", String(read(filepath, 3680))).captures[1]
sizey = match(r" SizeY=\"(.*?)\"", String(read(filepath, 3680))).captures[1]

println(framerate)
println(sizex)
println(sizey)