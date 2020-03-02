#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFileDialog>
#include <QDir>
#include "ExonChaining.h"
#include "LocalAlignment.h"
#include <stack>
#include <array>
#include <QDebug>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->tbl_result->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_btn_browseUnknownDNA_file_clicked()
{
    QString filters = "All Files (*.*) ;; Text Files (*.txt) ;; FASTA Files (*.fasta)";
    QString file_name = QFileDialog::getOpenFileName(this, "Open Unknown DNA file", QDir::currentPath(), filters);
    ui->txt_UnknownDNAS_path->setText(file_name);
}

void MainWindow::on_btn_browseDNATemplate_file_clicked()
{
    QString filters = "All Files (*.*) ;; Text Files (*.txt) ;; FASTA Files (*.fasta)";
    QString file_name = QFileDialog::getOpenFileName(this, "Open mRNA Template file", QDir::currentPath(), filters);
    ui->txt_mRNA_path->setText(file_name);
}

void MainWindow::on_btn_process_clicked()
{
    ui->tbl_result->clearContents();
    LocalAlignment la (ui->txt_UnknownDNAS_path->text().toStdString(),ui->txt_mRNA_path->text().toStdString());
    ExonChaining ec (la.getExons());
    stack<array<unsigned int, 3>> intervals = ec.get_intervals();
    int size = intervals.size();
    int i = 0;
    ui->tbl_result->setRowCount(size);
    while (!intervals.empty())
    {
        ui->tbl_result->setItem(i,0,  new QTableWidgetItem(QString::number(intervals.top()[0])));
        ui->tbl_result->setItem(i,1,  new QTableWidgetItem(QString::number(intervals.top()[1])));
        ui->tbl_result->setItem(i,2,  new QTableWidgetItem(QString::number(intervals.top()[2])));
        intervals.pop();
        i++;
    }
}
