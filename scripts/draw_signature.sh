#!/bin/bash
data_type=$1
number_of_topic=$2

if [ ${data_type} -eq 1 ]; then
    python Drawing/matching.py ${data_type} ${number_of_topic}
    python Drawing/draw_signature_pair.py ${data_type} ${number_of_topic}
    python Drawing/draw_arrangement.py ${data_type} ${number_of_topic}
fi

if [ ${data_type} -eq 2 ]; then
    python Drawing/convert_M2.py  ${number_of_topic}
    python Drawing/matching.py ${data_type} ${number_of_topic}
    python Drawing/draw_signature_pair.py ${data_type} ${number_of_topic}
    python Drawing/draw_arrangement.py ${data_type} ${number_of_topic}
    python Drawing/custom_draw_M2.py ${number_of_topic}
fi

