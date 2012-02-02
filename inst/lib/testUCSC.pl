#!/usr/bin/perl -w

use LWP;

##my $ARGV = @ARGV
my $ua = new LWP::UserAgent();
$ua->proxy(['http'] => "http://cache2.na.novartis.net:80");
my $resp = $ua->post("http://genome.ucsc.edu/cgi-bin/hgGateway");
my $code = $resp->code();
if ( !$resp->is_success ) {
   print STDERR "error: " . $resp->error_as_HTML . "\n";
}

my $hgsid;
if (scalar(@ARGV) < 4) {
    my $txt = $resp->content();
    ($hgsid) = $txt =~ /NAME\=\'hgsid\'\s+VALUE\=\'(\d+)\'/;
    print "$hgsid\n";
} else {
    $hgsid = $ARGV[3];
}
$req = $ua->post("http://genome.ucsc.edu/cgi-bin/hgCustom",
		 ["hgt.customFile" => [$ARGV[0]],
		  "org" => $ARGV[1],
		  "db" => $ARGV[2],
		  "hgsid"=>$hgsid,
		  "clade" => 'mammal'],
		 'Content_Type' => 'form-data');

$code = $req->code();
if ( !$req->is_success ) {
   print STDERR "error: " . $req->error_as_HTML . "\n";
}
