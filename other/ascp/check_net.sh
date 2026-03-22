sar -n DEV 1 | awk '/bond0/{ rx=$5*8/1000; tx=$6*8/1000; printf "RX: %.0f Mbps | TX: %.1f Mbps\n", rx, tx }'
