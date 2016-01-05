import sys
import os

header = ['Sample', 'Mouse_single', 'Mouse_multiple', 'Human', 'Other', 'Unmapped']
print '\t'.join(header)
for fi in sys.argv[1:]:
    sample = os.path.basename(fi).split('.')[0]
    if sample[-7:] == '_screen':
        sample = sample[:-7]
    with open(fi) as screen_results:
        results = {}
        for line in screen_results:
            fields = line.strip().split()
            if len(fields) and fields[0][0] != '#' and fields[0] != 'Library':
                if fields[0] == '%Hit_no_libraries:':
                    results['unmapped'] = int(float(fields[1]) / 100.0 * results['no_reads'])
                    continue
                results[fields[0] + '_single'] = int(fields[4])
                results[fields[0] + '_multiple'] = int(fields[6])
                results['no_reads'] = int(fields[1])
        
    if not len(results):
        data = ['0'] * 5

    else:

        try:
            data = [results['Mouse_single'],
                    results['Mouse_multiple'],
                    results['Human_single'] + results['Human_multiple']]
        except:
            sys.exit('Malformed file: {0}'.format(fi))

        data.append(results['no_reads'] - sum(data) - results['unmapped'])
        data.append(results['unmapped'])
        data = map(lambda i:str(i / float(sum(data))),data)


    data = [sample] + data
    print '\t'.join(data)
