
from subprocess import check_output as co
from datetime import datetime
from dateutil.relativedelta import relativedelta
import time

def add_months(sourcedate,months):
    month = sourcedate.month - 1 + months
    year = int(sourcedate.year + month / 12 )
    month = month % 12 + 1
    day = min(sourcedate.day,calendar.monthrange(year,month)[1])
    return datetime.date(year,month,day)

def f(sourcedate):
    return sourcedate.strftime('%Y-%m-%d')

def d(sourcedate):
    return time.mktime(sourcedate.timetuple())


start_date    = datetime(2014, 10, 1)
end_date      = datetime(2016, 8, 1)
current_date  = start_date
command       = 'git log --pretty=format:"%h%x09%an%x09%x09%s" --after="{range_start_s:s}" --before="{range_end_s:s}"'

while (d(current_date) < d(end_date)):
    range_start = current_date
    range_end = current_date + relativedelta(months=1)
    current_date = range_end
    range_start_s = f(range_start)
    range_end_s = f(range_end)
    print 'FROM {range_start_s:s} TO {range_end_s:s}'.format(**locals())
    cmd = command.format(**locals())
    lines = co(cmd.split())
    for line in lines.splitlines():
        if line.find('Jan_Hybs') == 9 or line.lower().find('x3mspeedy') == 9:
            print line
    


# print co('pwd')


# 
# 
# from dateutil.parser import parse
# dt = parse(start_date)
# # datetime.datetime(2010, 2, 15, 0, 0)
# print(dt.strftime('%Y-%m-%d'))

#$ git log --oneline --decorate --graph --after="2014-10-1" --before="2014-11-1"  --author="Jan_Hybs\|x3mSpeedy\|Jan"

