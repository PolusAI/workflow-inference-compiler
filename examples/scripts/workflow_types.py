# NOTE: Currently all types need to have a format, even if it is a placeholder.
string = {'type': 'string', 'format': 'edam:format_2330'}
integer = {'type': 'int', 'format': 'edam:format_2330'}
floating = {'type': 'float', 'format': 'edam:format_2330'}
# boolean not suported because CWL passes --flag, not --flag true
#boolean = {'type': 'boolean', 'format': 'edam:format_2330'}

textfile = {'type': 'File', 'format': 'edam:format_2330'} # textual format
csvfile = {'type': 'File', 'format': 'edam:format_3752'}
xlsxfile = {'type': 'File', 'format': 'edam:format_3620'}

pdbfile = {'type': 'File', 'format': 'edam:format_1476'}
mol2file = {'type': 'File', 'format': 'edam:format_3816'}
sdffile = {'type': 'File', 'format': 'edam:format_3814'}


pdbfiles = {'type': 'File[]', 'format': 'edam:format_1476'}
mol2files = {'type': 'File[]', 'format': 'edam:format_3816'}
sdffiles = {'type': 'File[]', 'format': 'edam:format_3814'}
