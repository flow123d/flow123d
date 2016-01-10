"""Static constants.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""

import platform

CONTEXT_NAME = 'ModelEditor'

DEFAULT_FONT = 'Courier,11,-1,5,50,0,0,0,0,0'
if platform.system() == 'Windows' and platform.release() != 'XP':
    DEFAULT_FONT = 'Consolas,11,-1,5,50,0,0,0,0,0'
