# The MIT License (MIT)
#
# Copyright (c) 2015 Kale Kundert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

white, black = '#ffffff', '#000000'
grey = '#eeeeec', '#d3d7cf', '#babdb6', '#888a85', '#555753', '#2e3436'

red =    '#ef2929', '#cc0000', '#a40000'
#red =    '', '#fa770b', '#cc0000', '#a40000'
orange = '#fcaf3e', '#f57900', '#ce5c00'
yellow = '#fce94f', '#edd400', '#c4a000'
green =  '#8ae234', '#73d216', '#4e9a06'
blue =   '#729fcf', '#3465a4', '#204a87'
#blue =   '','#00b0f0', '#3465a4', '#204a87'
purple = '#ad7fa8', '#75507b', '#5c3566'
brown =  '#e9b96e', '#c17d11', '#8f5902'

def from_cycle(index):
    cycle = (blue[1], red[1], green[2], orange[1], purple[1], brown[1],
             blue[0], red[0], green[1], orange[0], purple[0], brown[0])
    return cycle[index % len(cycle)]
