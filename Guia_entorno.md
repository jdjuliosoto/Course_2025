# Estadística multivariante en el análisis de datos genómicos 2025

## Entorno
Windows Subsystem for Linux versión 2 (WSL2) es una herramienta muy útil especialmente en el análisis genómico, ya que permite ejecutar herramientas bioinformáticas típicas de sistemas Unix/Linux directamente desde Windows.

*"permite a los desarrolladores instalar una distribución de Linux (como Ubuntu, OpenSUSE, Kali, Debian, Arch Linux, etc.) y usar aplicaciones, utilidades y herramientas de línea de comandos de Bash directamente en Windows, sin modificar, sin la sobrecarga de una máquina virtual tradicional o una configuración de arranque dual."*


Utilizar el subsistema de Windows para Linux tiene muchas ventajas para el análisis genómico. Pertime tener acceso nativo a herramientas bioinformáticas (sin necesidad de máquinas virtuales pesadas), lo que permite asignar mayor cantidad de recursos a los análisis.
Integración con Windows: puedes acceder a tus archivos locales desde Linux.
Uso de scripts bash, Python, R, etc., con compatibilidad total.
Entorno de desarrollo flexible y potente.

## 🐧 Instalar Ubuntu en Windows con WSL2
### Requisitos previos:
* Sistema operativo: Windows 10 o 11 (versión 64 bits)
* Tener permisos de administrador
* Espacio en disco disponible (recomendado al menos 10 GB)
* Habilitada la característica máquina virtual

### Paso 1: Habilitar WSL en su sistema
Abra PowerShell como administrador (menú > Inicio de PowerShell > , haga clic con el botón derecho en > Ejecutar como administrador) y escriba este comando:

```PowerShell
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
```
### Paso 2: Comprobación de los requisitos para ejecutar WSL 2
Para actualizar a WSL 2, debe ejecutar Windows 10 (revisar versiones en la página oficial de Microsoft) o Windows 11.
Para comprobar la versión y el número de compilación, seleccione la tecla del logotipo de Windows + R, escriba winver y seleccione OK. Actualice a la versión más reciente de Windows en el menú Configuración.

### Paso 3: Habilitación de la característica máquina virtual
Antes de instalar WSL 2, debe estar habilidatada la característica opcional Plataforma de máquina virtual. La máquina requerirá funcionalidades de virtualización para usar esta característica. Revisar (ctrl + shift + escape) 

En caso de no tener habilitada esta opción abra PowerShell como administrador y ejecute:

```PowerShell
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```
Reinicie la máquina para completar la instalación y actualización de WSL a WSL 2.

### Paso 4: Descarga del paquete de actualización del kernel de Linux

Descargue el paquete más reciente: https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi
Ejecute el paquete de actualización descargado en el paso anterior. (Haga doble clic para ejecutar; se le pedirán permisos elevados, seleccione "Sí" para aprobar esta instalación).

### Paso 5: Establecer WSL 2 como versión predeterminada
Abra PowerShell y ejecute este comando para establecer WSL 2 como versión predeterminada al instalar una nueva distribución de Linux:

```PowerShell
wsl --set-default-version 2
```
### Paso 6: Instalación de la distribución de Linux que prefiera
Abra Microsoft Store y seleccione Ubuntu 24.04.1 LTS. Descarguelo e instalelo. 
A continuación, deberá crear una cuenta de usuario y una contraseña para la nueva distribución de Linux.

Puede enumerar las distribuciones de Linux instaladas y comprobar la versión de WSL a la que está configurada cada una escribiendo el comando: wsl -l -v en PowerShell o en el símbolo del sistema de Windows.


Actualizar funciones basicas de Linux:

```

# Actualizar todos los paquetes
sudo apt update && sudo apt upgrade -y

# Incluye compiladores y herramientas necesarias para compilar codigo desde fuentes
sudo apt install build-essential

# Herramientas para descargar archivos desde la web.
sudo apt install curl wget

# Editores de texto para la linea de comandos.
sudo apt install vim nano

# Herramienta para monitorizar el uso de CPU, memoria y procesos en tiempo real.
sudo apt install htop
sudo apt install glances
```
Puede usar cualquier carpeta de Windows desde WSL accediendo a /mnt/c/ (por ejemplo):

```
cd /mnt/c/Users/usuario/Documents/proyecto_genomico
```

## Uso de repositorios para el control de versiones

En los últimos años, la forma en que trabajamos con código, datos y documentación ha evolucionado significativamente. En el ámbito del análisis de datos genómicos, donde los proyectos suelen involucrar múltiples archivos (scripts de R, archivos FASTQ, BAM, VCF, anotaciones, resultados intermedios, visualizaciones, etc.), es fundamental contar con una estructura organizada y reproducible.

Una herramienta clave para lograrlo es el uso de repositorios bajo control de versiones , y dentro de estas herramientas, Git se ha convertido en un estándar ampliamente adoptado. Combinado con plataformas como GitHub , Git permite no solo organizar mejor nuestros proyectos, sino también colaborar, compartir y mantener un historial claro de los cambios realizados a lo largo del tiempo.
¿Qué es un repositorio? Un repositorio (o repo , por sus siglas en inglés) es un espacio centralizado donde se almacenan los archivos de un proyecto junto con su historial completo de cambios .

Se puede pensar como una carpeta inteligente que no solo guarda los archivos actuales, sino también:

* Cada versión anterior de esos archivos,
* Quién hizo cada cambio,
* Cuándo se hizo,
* Y por qué se hizo (a través de mensajes de confirmación o commit messages ).

RStudio, por otro lado, facilita enormemente esta integración al permitirnos trabajar directamente con repositorios de GitHub desde su interfaz. Esto nos permite crear proyectos encapsulados que contienen todo lo necesario para ejecutar y replicar nuestro análisis: código, datos, documentación e incluso entornos personalizados.

Al vincular un proyecto de RStudio con un repositorio de GitHub, ganamos en varios aspectos fundamentales:

* Organización : Todo el proyecto vive en un único directorio, lo que facilita su manejo y portabilidad.
* Reproducibilidad : Garantizamos que cualquier persona pueda replicar el análisis exactamente como fue diseñado, gracias a que todos los archivos necesarios están disponibles y versionados.
* Control de cambios : Podemos rastrear cada modificación realizada, revertir errores y entender cómo ha evolucionado nuestro trabajo.
* Colaboración : Es fácil compartir el proyecto con otros investigadores o estudiantes, permitiendo comentarios, revisiones y contribuciones en tiempo real.
* Actualización continua : Si el repositorio pertenece a un curso, taller o biblioteca pública, podemos mantenernos al día incorporando nuevos contenidos o correcciones sin perder nuestro propio trabajo local.

A lo largo de este curso, se comentará cómo integrar estas prácticas en los proyectos de análisis genómico, comenzando por la creación de un proyecto en RStudio asociado a un repositorio de GitHub.

## Qué es Git?

Git es un software de código abierto para el control de versiones. Es una especie de administrador de proyectos no amigable.

### Instala Git
Una vez actualizado, instalar Git en Ubuntu con este comando:

```
sudo apt install git
```
Verificar que se instaló correctamente
```
git --version
# debería verse algo como: git version 2.xx.x
```

Configura Git (una sola vez): Es importante configurar el nombre de usuario y correo electrónico, ya que Git los usa para identificar quién hizo cada cambio.
```
git config --global user.name "SuNombre"
git config --global user.email "sucorreo@example.com"

# verificar configuración
git config --list
```

## Términos básicos 

**1. Repositorio (repo):**
Es la "carpeta inteligente" donde resides su proyecto.
Contiene sus archivos + el historial completo de cambios realizados.
Puede ser local (en tu computadora) o remoto (en GitHub, GitLab, etc.).

**2. Commit:**
Es una "instantánea" o versión guardada de su trabajo.
Cada commit tiene un mensaje que explica qué se hizo (git commit -m "mensaje").
Ejemplo: git commit -m "Añadido script de filtrado de SNPs"

**3. Branch (rama):**
Es una línea de desarrollo independiente.
La rama principal suele llamarse main o master. Se usan ramas para probar cosas nuevas sin afectar el proyecto principal.

Comandos útiles:
* git branch → ver ramas
* git checkout nombre_rama → cambiar de rama
* git switch -c nueva_rama → crear y cambiar a una nueva rama

**4. Merge:**
Es unir dos ramas. Por ejemplo, integrar cambios de una rama de prueba a la principal.
Ejemplo: después de trabajar en rama-analisis, puede hacer merge con main.

**5. Pull:**
Descargar cambios desde un repositorio remoto (como GitHub).
git pull origin main → trae las últimas actualizaciones del repositorio remoto.

**6. Push:**
Subir sus commits locales al repositorio remoto.
git push origin main → envía sus cambios a GitHub.

**7. Clone:**
Copiar un repositorio remoto a su computadora
git clone URL_DEL_REPO → crea una copia local del proyecto.

**8. Status:**
Muestra el estado actual del repositorio: qué archivos han cambiado, si hay commits pendientes, etc.
git status → es muy útil antes de hacer un commit.

**9. Add:**
Añade cambios al área de preparación (staging ) antes de hacer commit.
git add nombre_archivo.R → prepara solo ese archivo.
git add . → prepara todos los archivos modificados.

**10. Log:**
Muestra el historial de commits realizados.
git log → le muestra quién hizo qué y cuándo.

### Resumen visual

| Acción         | Comando típico                      | Descripción                                      |
|----------------|-----------------------------------|--------------------------------------------------|
| Clonar         | `git clone URL`                   | Copiar un repositorio remoto a su máquina        |
| Añadir         | `git add .`                       | Preparar todos los cambios para hacer commit     |
| Guardar        | `git commit -m "mensaje"`         | Crear una nueva versión del proyecto             |
| Ver estado     | `git status`                      | Ver qué archivos han sido modificados            |
| Subir cambios  | `git push origin main`            | Enviar sus commits al repositorio remoto         |
| Descargar      | `git pull origin main`            | Obtener los últimos cambios del repositorio remoto |
| Ver historial  | `git log`                         | Revisar el historial de cambios                  |
| Ver ramas      | `git branch`                      | Mostrar las ramas existentes                     |
| Cambiar rama   | `git checkout nombre_rama`        | Moverse entre ramas                              |
| Crear rama     | `git switch -c nueva_rama`        | Crear y moverse a una nueva rama                 |
| Elimina rama   | `git branch -d nombre_rama`       | Elimina rama                                     |


### Commits y Branches en Git
A, B, C: Commits realizados en la rama main. En el commit C, se crea una nueva rama llamada feature-analisis. D, E: Nuevos commits realizados únicamente en la rama feature-analisis.
```
main (rama principal)
|
A---B---C
         \
          \___ feature-analisis (rama secundaria)
               \
                D---E
```
Se puede hacer un merge de feature-analisis con main para integrar los cambios, o eliminar la rama si no se usará
```
main (rama principal)
|
A---B---C---------------------------F (merge)
         \                         /
          \___ feature-analisis __/
               \
                D---E
```
Aquí el commit F representa el resultado del merge entre main y feature-analisis.

### Las 3 etapas de un archivo en Git

Git maneja un sistema de tres estados principales para los archivos dentro de un repositorio:
* Modified (Modificado): El archivo ha sido editado, pero Git aún no lo conoce.
* Staged (Preparado): El archivo se ha añadido al "stage" (área de preparación) para incluirse en el próximo commit.
* Committed (Confirmado): El archivo ha sido guardado permanentemente en la base de datos de Git.

```
# 1. Edita un archivo (por ejemplo script.R)
nano script.R

# 2. Verifica el estado del repositorio
git status

# 3. Añade los cambios al staging area
git add script.R

# 4. Confirma los cambios (commit)
git commit -m "Añadida función de filtrado de SNPs"
```

```
┌──────────────┐     git add      ┌───────────────┐     git commit     ┌───────────────┐
│              │  ────────────>   │               │  ─────────────>    │               │
│  Modified    │     (staging)    │    Staged     │     (commited)     │   Committed   │
│              │  <───────────    │               │                    │               │
└──────────────┘     git reset    └───────────────┘                    └───────────────┘
```

## Git en windows
Git para Windows, también conocido como msysgit o "Git Bash", permite obtener Git además de otras herramientas útiles, como la shell Bash. 
Git para Windows deja el ejecutable de Git en una ubicación convencional, lo que le ayudará y a otros programas, como RStudio, a encontrarlo y usarlo. 
Esto también facilita la transición a un uso más experto, ya que la shell "Git Bash" le resultará útil al explorar fuera de R/RStudio.

Para instalarlo, descargar el ejecutable de https://git-scm.com/downloads e instalarlo con la configuración predeterminada.

### Instalar R y rstudio
R es un entorno de programación diseñado específicamente para el análisis estadístico y la visualización de datos, lo que lo convierte en una herramienta ideal para científicos, investigadores y analistas que necesitan procesar grandes volúmenes de información con precisión. 
Su amplia comunidad y su ecosistema de paquetes permiten acceder a métodos avanzados de manera sencilla y reproducible. RStudio , por su parte, ofrece una interfaz intuitiva y completa que facilita el trabajo con R, integrando herramientas clave como control de versiones, edición de scripts, visualización interactiva y generación de informes dinámicos. 
Juntos, forman una plataforma poderosa, flexible y accesible, perfecta para explorar, analizar y comunicar resultados a partir de datos complejos, como los que se manejan en estudios genómicos y bioinformáticos.

Para instalar R, descargar el ejecutable correspondiente de https://cran.r-project.org e instalarlo en windows ya que se utilizará en Windows.
Para insalar rstudio descargar el ejecutable de https://posit.co/downloads/ e instalarlo en Windows.

RStudio viene con una versión portable de Git para facilitar el uso del control de versiones directamente desde la interfaz. Para utilizarla, abrir el ejecutable de rstudio. 
Una vez abierto dirigirse al panel de la consola (inferior izquierdo) y acceder a la pestaña terminal. En la terminal acceder a la version instalada de git con el comando: 

```
which git
git --version
```
En este caso también hay que configurar Git con un nombre de usuario y un correo electrónico. Para hacer esto moverse a la pestaña consola y escribir en esta el siguiente comando: 

```
install.package("usethis")
library(usethis)
usethis::edit_git_config()
# Modificar en el fichero ".gitconfig" los apartados: "name" y "email" 
# y guardar el fichero y puede cerrarse
```

### Inicializando un nuevo repositorio en Rstudio

Una de las grandes ventajas de trabajar con RStudio es la posibilidad de crear proyectos asociados directamente a un repositorio de GitHub, lo cual se logra mediante la clonación del repositorio dentro de un nuevo proyecto en RStudio. 
Esta integración permite aprovechar al máximo las funcionalidades de ambos entornos: por un lado, el proyecto queda perfectamente organizado dentro de RStudio, encapsulando todos los archivos relevantes (scripts, datos, resultados, documentación) en una misma unidad de trabajo; 
por otro lado, se encuentra bajo control de versiones gracias a Git y alojado en GitHub, lo cual garantiza trazabilidad, seguridad y capacidad de colaboración. 
Esto resulta especialmente valioso cuando trabajamos con repositorios mantenidos por terceros, ya que nos permite mantener nuestro código actualizado fácilmente incorporando los últimos cambios desde el repositorio original. 
Además, si contamos con permisos de escritura, podemos realizar modificaciones locales y sincronizarlas directamente hacia el repositorio remoto desde la propia interfaz de RStudio, facilitando un flujo de trabajo ágil, transparente y altamente colaborativo. 

La función use_git() agregará un repositorio Git (a menudo denominado “repositorio”) a un proyecto RStudio existente.

En RStudio:
1. Crear un nueva proyecto
2. Seleccionar “Nuevo Directorio”
3. Proyecto
4. Activar: “Create a git repository”

Otra manera de crear un repositorio es mediante la librería *usethis* por medio de la consola.
```
library(usethis)
usethis::use_git()
# Elegir siempre la opción: 1
# Y ante la ventana, seleccionar: "Save"
```
Si todo salio bien, una pestaña con el nombre de Git aparecerá en la ventana superior derecha. Al seleccionarla se puede visitar todo el historial del proyecto al dar click en el ícono del reloj.
Cada vez que se desee realizar un *commit* se debe dar click en el ícono con el nombre commit. Se selecciona el scritp que quiere confirmarse, se agrega el descriptor del cambio y se confirma.
De esta manera se puede revisitar cada cambio realizado durante el proceso de edición.

---
Pueden ampliar esta breve explicación mediante alguno de los muchos materiales que pueden encontrar por internet.
* https://learn.microsoft.com/es-es/windows/wsl/install
* https://swcarpentry.github.io/git-novice-es/
* https://git-scm.com/book/en/v2
* https://cran.r-project.org

