#include <stdio.h>
#include <iostream>
#include <gtkmm.h>
#include <chrono>
#include <cmath>
#include "headers/worker.hpp"

Gtk::ProgressBar *progress_bar;
Gtk::Label* eta;
std::thread* worker_thread;
Glib::RefPtr<Gtk::Application> app;

void notify(float fraction, std::chrono::milliseconds time, int runs_left) {
	progress_bar->set_fraction(fraction);
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(time * runs_left);
	int count = duration.count();
	int minutes = 0;
	int seconds = 0;
	if(count > 60) {
		minutes = floor(count / 60);
		seconds = count % 60;
	}
	if(minutes == 0) {

		eta->set_text("Estimated time left: " + std::to_string(seconds) + " seconds");
	} else {
		eta->set_text("Estimated time left: " + std::to_string(minutes) + " minutes, " + std::to_string(seconds) + " seconds");
	}
	if(fraction == 1.0) {
		app->quit();
	}
}

void worker_code() {
	worker new_worker;//What a load of sorcery.
	new_worker.do_work(notify);
}

int main(int argc, char** argv) {
	app = Gtk::Application::create(argc,argv, "org.gtkmm.lbm_simulation");
	Gtk::Window* window = nullptr;
	Glib::RefPtr<Gtk::Builder> builder = Gtk::Builder::create_from_file("gui.glade");
	builder->get_widget("mainWindow", window);
	progress_bar = nullptr;
	eta = nullptr;
	builder->get_widget("progressBar", progress_bar);
	builder->get_widget("eta", eta);
	window->set_title("LBM simulation");
  window->set_keep_above(true);
	worker_thread = new std::thread(worker_code);
	return app->run(*window);
}
