import gamtools

def test_members():

    for module in ('call_windows', 'compaction', 'cosegregation', 'enrichment', 'matrix',
                   'permutation', 'resolution', 'radial_position', 'segregation' ,'utils'):
        print(dir(getattr(gamtools, module)))
