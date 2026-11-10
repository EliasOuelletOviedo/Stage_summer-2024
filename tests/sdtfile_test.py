from sdtfile import SdtFile

filepath = "/Users/eliasouellet-oviedo/Documents/Stages/Stage2/Partie2/Données/2026-03-12/FOVs/fov_3.sdt"

s = SdtFile(filepath)

print(s.header.revision)
print(s.block_headers[0])
print(s.block_headers[0].block_type)
print(type(s.data[0]))
print(s.data[0].dtype)
print(s.data[0].shape)