import sys
results = {}
with open(sys.argv[1]) as screen_results:
    for line in screen_results:
        fields = line.strip().split()
        if len(fields) and fields[0][0] != '#' and fields[0] != 'Library':
            if fields[0] == '%Hit_no_libraries:':
                results['unmapped'] = int(float(fields[1]) / 100.0 * results['no_reads'])
                continue
            results[fields[0] + '_single'] = int(fields[4])
            results[fields[0] + '_multiple'] = int(fields[6])
            results['no_reads'] = int(fields[1])
header = ['Mouse_single', 'Mouse_multiple', 'Human', 'Other', 'Unmapped']
data = [results['Mouse_single'],
        results['Mouse_multiple'],
        results['Human_single'] + results['Human_multiple']]
data.append(results['no_reads'] - sum(data) - results['unmapped'])
data.append(results['unmapped'])
print '\t'.join(header)
print '\t'.join(map(lambda i:str(i / float(sum(data))),data))
