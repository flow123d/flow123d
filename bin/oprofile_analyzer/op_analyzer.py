"""
Selective profiling via. statistical OProfile profiler and manual blowup (repetition) of selected
part of code.

Usage: op_analyzer.py  profiling_dir

Hardwired is also the base_dir.
Both folders contain two files mass_1.txt nad mass_500.txt.
It is assumed that the first is a reference OProfile report while the second have 500x repetition of a sleceted part of code.

Output:

     0    0.048s   2.615s |  0.013%   2.364% | TOTAL [s]
  8.45%    8.45%    4.16% |   0.00%    1.55% |  dgemv_
 15.77%    7.32%    0.23% |   0.00%    0.00% |  MassAssemblyDG<2u, ConcentrationTransportModel>::cell_integral(unsigned int, unsigned int)
 20.59%    4.82%    3.48% |   0.00%    0.00% |  MapPiola<3u>::update_values(FEValues<3u>&, ElementValues<3u> const&, FEValues<3u>::FEInternalData const&)
 24.84%    4.25%    1.56% |   0.00%    1.16% |  ilaenv_

1. col - cumulative percentage of the 2. col
2. col - percentage of time consumed by a symbol within the blown part ( time in symbol in blown / time of blown)
3. col - percentage of time consumed by a symbol out of the blown part ( time in symbol out of blown / time out of blown)
4. and 5. - same as 2. and 3. but for the base directory (assumed a base or reference commit)
6. col the symbol name

High percentage in 2. col and small in 4. col hints a performance regression in the blown code part. Unnecessary call of the symbol.
"""
import re
import os
import sys



def symbol_pair(line):

    cols = line[:-1].split('  ')
    #print(cols)
    samples, percent, line_nr, image, app, symbol = [t for t in cols if t]
    factor = float(percent) / 100
    return symbol, factor


def read_oprofile(f):
    with open(f, "r") as f:
        line = f.readline()
        #print(line)
        m = re.match("TEST TIME: ([0-9.]*) s", line)
        total_time = float(m.groups()[0])
        #print(total_time)
        
        line = f.readline()
        #print(line)
        #assert(line == "\n")
        
        #header
        line = f.readline()
        
        symbols = [symbol_pair(line) for line in f if line!='\n']
        return total_time, dict(symbols)

def process_dir(dir):
    blowup = 500
    total_time_r, symbols_r = read_oprofile(os.path.join(dir, "mass_1.txt"))
    total_time_b, symbols_b = read_oprofile(os.path.join(dir, "mass_500.txt"))
    symbols = {}
    sum_frac_r = 0
    sum_frac_b = 0

    def real_times(time_frac_b, time_frac_r):
        """
        Input: relative total time consumed by symbol in: blown variant, reference variant
        Return absolute times consumed by the symbol in: blown part, remaining reference part.
        """
        if time_frac_b * total_time_b > time_frac_r * total_time_r:
            blow_time = (time_frac_b * total_time_b - time_frac_r * total_time_r) / (blowup - 1)
            ref_time = (blowup * time_frac_r * total_time_r - time_frac_b * total_time_b) / (blowup - 1)
        else:
            blow_time = 0
            ref_time = time_frac_b * total_time_b
        return blow_time, ref_time

    for symbol, time_frac_b in symbols_b.items():
        try:
            time_frac_r = symbols_r.get(symbol, 0.0)
        except KeyError:
            continue


        sum_frac_b += time_frac_b
        sum_frac_r += time_frac_r
        symbols[symbol] = real_times(time_frac_b, time_frac_r)

    # compute REMINDER
    reminder_b = 1 - sum_frac_b
    reminder_r = 1 - sum_frac_r
    #blow_time = (reminder_b * total_time_b - reminder_r * total_time_r) / (blowup - 1)
    #ref_time = (blowup * reminder_r * total_time_r - reminder_b * total_time_b) / (blowup - 1)
    symbols['REMINDER'] = real_times(reminder_b, reminder_r)

    # total time
    # blow_time = (total_time_b - total_time_r) / (blowup - 1)
    # ref_time = (blowup * total_time_r - total_time_b) / (blowup - 1)
    symbols['TOTAL'] = real_times(1, 1)

    return  symbols




def format_table(base, new):

    def get_row(key):
        new_row = new[key]
        base_row = base.get(key, (0, 0))
        return new_row[0], new_row[1], base_row[0], base_row[1]

    
    totals = get_row('TOTAL')
    #print(totals)
    total_tb_new, total_tr_new, total_tb_base, total_tr_base = totals

   

    def format_row(tcumul, tb_new, tr_new, tb_base, tr_base, symbol):
        line = (f"{tcumul / total_tb_new * 100:6.2f}% | "
               f"{tb_new / total_tb_new * 100:6.2f}%  "
               f"{tr_new / total_tr_new * 100:6.2f}% | "
               f"{tb_base / total_tb_base * 100:6.2f}%  "
               f"{tr_base / total_tr_base * 100:6.2f}% | "
               f"{symbol}")
        print(line)


    rows = [(*get_row(k), k) for k in new.keys()]
    rows.sort(reverse = True)
    print( "       | time new: TB  TR | time base: TB TR | ")
    print(
        f'{0:6} | {total_tb_new:6.3f}s  {total_tr_new:6.3f}s | {total_tb_base:6.3f}s  {total_tr_base:6.3f}s | TOTAL [s]')
    cumul_blow_new = 0
    print( "=================================================== ")
    print( " cum. % | % new: TB     TR | % base: TB    TR | symbol")

    assert rows[0][4] == 'TOTAL', rows[0]
    for row in rows[1:]:
        cumul_blow_new += row[0]
        format_row(cumul_blow_new, *row)



if __name__ == "__main__":
    new_dir = sys.argv[1]
    base_dir = "6d25_304"
    #new_dir = "3b79"

    base_symbols = process_dir(base_dir)
    new_symbols = process_dir(new_dir)

    format_table(base_symbols, new_symbols)

