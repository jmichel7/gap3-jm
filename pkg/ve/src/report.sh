makemess() {
    echo Subject: NME build report
    echo `date` at `hostname` by $USER
}

makemess | /usr/bin/mail sal@cs.st-and.ac.uk

