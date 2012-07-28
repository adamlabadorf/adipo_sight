#!/usr/bin/python

import cgi
import cgitb
cgitb.enable()

import logging
import re


log_fn = 'adipo_site.log'
log_fmt = '%(asctime)s|[%(levelname)s] - %(message)s'
logging.basicConfig(format=log_fmt,filename=log_fn,level=logging.DEBUG)
debug = logging.debug
info = logging.info
warn = logging.warning
error = logging.error

from jinja2 import Environment, PackageLoader, Template

def log_to_html(l) :
    try :
        m = re.search('^(?P<time>.*)\|\[(?P<level>.*)\] - (?P<message>.*)$',l)# - (?P<message>.*)$',l)
        l = '<span class="%s">%s</span>'%(m.group('level').lower(),l)
    except :
        l = '<span class="error">Error parsing log line: %s</span>'%l
    return l

def get_log_content() :
    loader = PackageLoader('adipo_sight','data/tmpl')
    env = Environment(loader=loader)
    template = env.get_template('log.html')

    env.globals['log_to_html'] = log_to_html

    log = open(log_fn).readlines()[::-1]
    return template.render(log=log[:100])

if __name__ == '__main__' :

    print "Content-Type: text/html"
    print

    print get_log_content()

