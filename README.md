# Projekt

### plot_framework
For now (Mar 3rd) it includes just a simple `plot_attributes`-function. You can run it after plotting anything to customize the attributes of the plot.
Include it by 
```jl
include("plot_framework")
```


## Allgemeine Dinge zu `git`

* ES WIRD **NIEMALS** IN DEN UPSTREAM MASTER GEPUSHT!!!! ES WIRD GENAU GENOMMEN SOGAR **NIEMALS** IN IRGENDEINEN UPSTREAM BRANCH GEPUSHT, WENN ES MEHR ALS NUR DEN UPSTREAM MASTER GIBT!!!!

### Erstes setup:

1.  Dieses Repository in den eigenen Namespace forken (siehe `fork button`)
2.  Auf lokalem Rechner folgende Befehle ausführen (der angegebene upstream Link ist jetzt der für `ssh`-Login!):
```sh
git  clone <url_von_eurer_fork>
git remote add upstream git@gitlab.gwdg.de:data_assimilation_group/data_assimilation.git 
git fetch --all
```
### Workflow:
* Vom `upstream` aktuelle Version auf lokalen Rechner pullen 
```sh 
git pull upstream <branch>
```
* Code lokal weiterentwickeln
* README lokal updaten: kurze Zusammenfassung des features, das ihr neu eingebaut habt.
* Check: Seid ihr lokal auf dem richtigen branch? Mit diesem Befehl bekommt ihr angezeigt, auf welchem branch ihr euch lokal befindet, welche Dateien verändert wurden und welche evtl noch gar nicht von `git` erfasst wurden.
```sh
git status
```
* Veränderte Dateien für den commit vorbereiten
```sh
git add <dateiname>
```
* `commit`: Eine Art Zwischenspeicherung. Commits immer mit kurzer Beschreibung, was ihr verändert habt, ausstatten:
```sh
git commit -m"<commit_message>"
```
* Wenn ihr mit eurer Arbeit soweit zufrieden seid: in `origin` pushen (Origin ist eure eigene Fork). Das macht ihr **nicht jedes Mal, wenn ihr committed**.
```sh
git push origin <branch>
```
* Merge request mit `upstream` erstellen (online, also auf der website). Nur ich (Inga) als Besitzerin dieses Repositorys habe die Möglichkeit, diese merge requests zu bewilligen.



### Weitere nützliche Befehle
* Ihr könnt sehen, welche branches es gibt, indem ihr folgenden Befehl eingebt:
```sh
git branch
```
* Einen neuen branch erstellt ihr mit
```sh
git checkout -b <neuer_branch_name>
```
* Falls ihr nicht auf dem richtigen branch seid: in den branch gehen, in den ihr eure aktuelle Arbeit hochladen wollt:
```sh
git checkout <branch_name>
```
* Damit unnötige files (zum Beispiel ausführbare Dateien, Dateien mit Daten zum Plotten, etc) **nicht** das Repository zumüllen, könnt ihr sie zu `.gitignore` hinzufügen. Jede Zeile in dieser Datei gibt eine Datei(enart) an, die nicht von `git` getrackt werden soll:
```sh
cat .gitignore
# ignore generated html files,
*.html
# except foo.html which is maintained by hand
!foo.html
```

### Einen `ssh`-key einrichten
* Zunächst einen key generieren mit
```sh
ssh-keygen -t ed25519 -C "vorname.nachname@stud.uni-goettingen.de"
```
* Den key ausgeben lassen mit
```sh
cat /home/nutzername/.ssh/id_ed25519.pub
```
* Die Ausgabe auf dem Terminal in die Zwischenablage kopieren 
* Online auf www.gitlab.gwdg.de oben ganz rechts unter usr -> Settings gehen, dann in der linken Leiste auf "SSH-keys".
* Den key in das Feld aus der Zwischenablage einfügen
* Gegebenenfalls einen Titel angeben
* Auf `add key` gehen.
* Fertig!

Mit diesem `ssh`-key könnt ihr nun euer `git`-Repository runterladen. Dabei fürs clonen einfach den oberen Link (ssh) kopieren und ganz normal im Terminal
```sh
git clone <url>
```
eingeben.

