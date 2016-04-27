"""Analysis related utilities.

.. codeauthor:: Tomas Krizek <tomas.krizek1@tul.cz>
"""


def replace_params_in_file(src_path, dst_path, params):
    with open(src_path, 'r') as src:
        text = src.read()
    text = replace_params_in_text(text, params)
    with open(dst_path, 'w') as dst:
        dst.write(text)


def replace_params_in_text(text, params):
    """
    Replace all parameters in input text with provided values.

    :param text: input text
    :param params: dictionary of param: value
    :return: text with replacements
    """
    for name, value in params.items():
        text = text.replace('<%s>' % name, str(value))
    return text

