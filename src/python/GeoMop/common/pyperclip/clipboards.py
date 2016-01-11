import sys
import subprocess

PY2 = sys.version_info[0] == 2


def init_xclip_clipboard():
    def copy_xclip(text):
        p = subprocess.Popen(['xclip', '-selection', 'c'], stdin=subprocess.PIPE, close_fds=True)
        p.communicate(input=text.encode('utf-8'))

    def paste_xclip():
        p = subprocess.Popen(['xclip', '-selection', 'c', '-o'], stdout=subprocess.PIPE, close_fds=True)
        stdout, stderr = p.communicate()
        return stdout.decode('utf-8')

    return copy_xclip, paste_xclip


def init_no_clipboard():
    class ClipboardUnavailable(object):
        def __call__(self, *args, **kwargs):
            pass

        if PY2:
            def __nonzero__(self):
                return False
        else:
            def __bool__(self):
                return False

    return ClipboardUnavailable(), ClipboardUnavailable()
